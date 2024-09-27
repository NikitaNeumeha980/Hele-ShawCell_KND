#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>

using namespace std;

int main(int argc,char* argv[])
{
    cout << "Permeability tensor calculation..." << endl;

    ofstream
        permeabilityTensorI1I2I3(static_cast<string>(argv[1]) + "/permeabilityTensorI1I2I3");

    ifstream
        mAveragedGradPDivNu(static_cast<string>(argv[1]) + "/mAveragedGradPDivNu"),
        averagedU(static_cast<string>(argv[1]) + "/averagedU");

    string
        string_buffer;

    int
        iteration_index(0),
        N(0);

    double
        *x,
        *x_prev,
        *rhs_buffer,
        *rhs,
        *matrix_components,
        **matrix_buffer,
        **matrix,
        diagonal_component,
        sum(0),
        sum_err(0),
        eps(1e-12);

    while
    (
        getline(mAveragedGradPDivNu, string_buffer)
    )
        N++;

    x = new double [N];
    x_prev = new double [N];
    rhs_buffer = new double [N];
    rhs = new double [N];
    matrix_components = new double [N];
    matrix_buffer = new double* [N];
    matrix = new double* [6];
    for
    (
        int i(0);
        i < N;
        i++
    )
        matrix_buffer[i] = new double [6];

    for
    (
        int i(0);
        i < 6;
        i++
    )
        matrix[i] = new double [6];


    if (!mAveragedGradPDivNu.is_open())
        cout << "File can not be open!" << endl;
    else
    {
        mAveragedGradPDivNu.clear();
        mAveragedGradPDivNu.seekg(0);

        averagedU.clear();
        averagedU.seekg(0);

        for
        (
            int i(0);
            i < N;
            i++
        )
        {
            getline(mAveragedGradPDivNu, string_buffer);
            matrix_components[i] = strtod(string_buffer.c_str(), 0);

            getline(averagedU, string_buffer);
            rhs_buffer[i] = strtod(string_buffer.c_str(), 0);

            //cout << matrix_components[i] << "\t" << rhs_buffer[i] << endl;
        }

        mAveragedGradPDivNu.close();
        averagedU.close();
    }

    for
    (
        int i(0);
        i < 6;
        i++
    )
    {
        x[i] = 0;
        x_prev[i] = 0;
        for
        (
            int j(0);
            j < 6;
            j++
        )
            matrix[i][j] = 0;
    }

    for
    (
        int i(0);
        i < N;
        i++
    )
    {
        for
        (
            int j(0);
            j < 6;
            j++
        )
            matrix_buffer[i][j] = 0;
    }

    matrix_buffer[0][0] = matrix_components[0];
    matrix_buffer[0][1] = matrix_components[1];
    matrix_buffer[0][2] = matrix_components[2];

    matrix_buffer[1][1] = matrix_components[0];
    matrix_buffer[1][3] = matrix_components[1];
    matrix_buffer[1][4] = matrix_components[2];

    matrix_buffer[2][2] = matrix_components[0];
    matrix_buffer[2][4] = matrix_components[1];
    matrix_buffer[2][5] = matrix_components[2];

    matrix_buffer[3][0] = matrix_components[3];
    matrix_buffer[3][1] = matrix_components[4];
    matrix_buffer[3][2] = matrix_components[5];

    matrix_buffer[4][1] = matrix_components[3];
    matrix_buffer[4][3] = matrix_components[4];
    matrix_buffer[4][4] = matrix_components[5];

    matrix_buffer[5][2] = matrix_components[3];
    matrix_buffer[5][4] = matrix_components[4];
    matrix_buffer[5][5] = matrix_components[5];

    matrix_buffer[6][0] = matrix_components[6];
    matrix_buffer[6][1] = matrix_components[7];
    matrix_buffer[6][2] = matrix_components[8];

    matrix_buffer[7][1] = matrix_components[6];
    matrix_buffer[7][3] = matrix_components[7];
    matrix_buffer[7][4] = matrix_components[8];

    matrix_buffer[8][2] = matrix_components[6];
    matrix_buffer[8][4] = matrix_components[7];
    matrix_buffer[8][5] = matrix_components[8];

    for
    (
        int i(0);
        i < 6;
        i++
    )
    {
        sum = 0;

        for
        (
            int j(0);
            j < N;
            j++
        )
            sum = sum + matrix_buffer[j][i] * rhs_buffer[j];

        rhs[i] = sum;

        for
        (
            int k(0);
            k < 6;
            k++
        )
        {
            sum = 0;

            for
            (
                int j(0);
                j < N;
                j++
            )
                sum = sum + matrix_buffer[j][i] * matrix_buffer[j][k];

            matrix[i][k] = sum;
        }
    }

    for
    (
        int i(0);
        i < 6;
        i++
    )
    {
        rhs[i] = rhs[i]/matrix[i][i];
        diagonal_component = matrix[i][i];

        for
        (
            int j(0);
            j < 6;
            j++
        )
            if(j == i)
                matrix[i][j] = 0;
            else
                matrix[i][j] = - matrix[i][j]/diagonal_component;
    }

    iteration_index = 0;

    do
    {
        for
        (
            int i(0);
            i < 6;
            i++
        )
        {
            sum = 0;
            x_prev[i] = x[i];

            for
            (
                int j(0);
                j < 6;
                j++
            )
                sum = sum + matrix[i][j]*x[j];

            x[i] = rhs[i] + sum;
        }

        sum = 0;
        sum_err = 0;

        for
        (
            int i(0);
            i < 6;
            i++
        )
        {
            sum = sum + x[i] * x[i];
            sum_err = sum_err + (x[i] - x_prev[i]) * (x[i] - x_prev[i]);
        }

        iteration_index += 1;

        cout << "Iteartion = " << iteration_index << ", eps = " << sqrt(sum_err / sum) << endl;

    }while
    (
        (sqrt(sum_err / sum) > eps)
    );

    permeabilityTensorI1I2I3
        << x[0] << " " << x[1] << " " << x[2] << endl
        << x[1] << " " << x[3] << " " << x[4] << endl
        << x[2] << " " << x[4] << " " << x[5] << endl
        << x[0] + x[3] + x[5] << endl
        << -(x[3] + x[0]) * x[5] + x[4] * x[4] - x[0] * x[3] + x[2] * x[2] + x[1] * x[1] << endl
        << x[0] * (x[3] * x[5] - x[4] * x[4]) - x[1] * (x[1] * x[5] - x[2] * x[4]) + x[2] * (x[1] * x[4] - x[2] * x[3]) << endl;

    delete [] x;
    delete [] x_prev;
    delete [] rhs_buffer;
    delete [] rhs;
    delete [] matrix_components;
    for
    (
        int i(0);
        i < N;
        i++
    )
        delete [] matrix_buffer[i];
    delete matrix_buffer;
    for
    (
        int i(0);
        i < 6;
        i++
    )
        delete [] matrix[i];
    delete matrix;
    
    cout << "ok" << endl;
    
    return 0;
}
 
