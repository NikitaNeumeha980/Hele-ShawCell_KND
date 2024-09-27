# state file generated using paraview version 5.4.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2229, 996]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0, 0.0, -4.999999873689376e-06]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.0, 0.0, 0.016387316371567536]
renderView1.CameraFocalPoint = [0.0, 0.0, -4.999999873689376e-06]
renderView1.CameraParallelScale = 0.00289778271311716
renderView1.Background = [0.32, 0.34, 0.43]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'OpenFOAMReader'
case1foam = OpenFOAMReader(FileName='/media/alexshtil/STORAGE/calculations/OpenFOAM/regularPorousMedia/twoPhase/incompressible/lattice/case2/case2.foam')
case1foam.SkipZeroTime = 0
case1foam.CaseType = 'Decomposed Case'
case1foam.MeshRegions = ['internalMesh']

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.05, 0.865003, 0.865003, 0.865003, 0.1, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'alphaphase1'
alphaphase1LUT = GetColorTransferFunction('alphaphase1')
alphaphase1LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
alphaphase1LUT.ColorSpace = 'RGB'
alphaphase1LUT.NanColor = [1.0, 0.0, 0.0]
alphaphase1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'alphaphase1'
alphaphase1PWF = GetOpacityTransferFunction('alphaphase1')
alphaphase1PWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'p_rgh'
p_rghLUT = GetColorTransferFunction('p_rgh')
p_rghLUT.RGBPoints = [-100.0, 0.231373, 0.298039, 0.752941, 200.0, 0.865003, 0.865003, 0.865003, 500.0, 0.705882, 0.0156863, 0.14902]
p_rghLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'p_rgh'
p_rghPWF = GetOpacityTransferFunction('p_rgh')
p_rghPWF.Points = [-100.0, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]
p_rghPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from case1foam
case1foamDisplay = Show(case1foam, renderView1)
# trace defaults for the display properties.
case1foamDisplay.Representation = 'Surface'
case1foamDisplay.ColorArrayName = ['POINTS', 'p_rgh']
case1foamDisplay.LookupTable = p_rghLUT
case1foamDisplay.OSPRayScaleArray = 'U'
case1foamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
case1foamDisplay.SelectOrientationVectors = 'U'
case1foamDisplay.ScaleFactor = 0.0006000000052154065
case1foamDisplay.SelectScaleArray = 'None'
case1foamDisplay.GlyphType = 'Arrow'
case1foamDisplay.GlyphTableIndexArray = 'None'
case1foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
case1foamDisplay.PolarAxes = 'PolarAxesRepresentation'
case1foamDisplay.ScalarOpacityFunction = p_rghPWF
case1foamDisplay.ScalarOpacityUnitDistance = 0.0001460986515426704
case1foamDisplay.GaussianRadius = 0.00030000000260770325
case1foamDisplay.SetScaleArray = ['POINTS', 'alpha.phase1']
case1foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
case1foamDisplay.OpacityArray = ['POINTS', 'alpha.phase1']
case1foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
case1foamDisplay.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for alphaphase1LUT in view renderView1
alphaphase1LUTColorBar = GetScalarBar(alphaphase1LUT, renderView1)
alphaphase1LUTColorBar.Title = 'marker'
alphaphase1LUTColorBar.ComponentTitle = ''
alphaphase1LUTColorBar.TitleFontSize = 25
alphaphase1LUTColorBar.LabelFontSize = 25
alphaphase1LUTColorBar.AutomaticLabelFormat = 0
alphaphase1LUTColorBar.LabelFormat = '%-#6.1f'
alphaphase1LUTColorBar.RangeLabelFormat = '%-#6.1f'
alphaphase1LUTColorBar.ScalarBarThickness = 50
alphaphase1LUTColorBar.ScalarBarLength = 0.3

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Title = '|U|, [m/s]'
uLUTColorBar.ComponentTitle = ''
uLUTColorBar.TitleFontSize = 25
uLUTColorBar.LabelFontSize = 25
uLUTColorBar.AutomaticLabelFormat = 0
uLUTColorBar.LabelFormat = '%-#6.2f'
uLUTColorBar.RangeLabelFormat = '%-#6.2f'
uLUTColorBar.ScalarBarThickness = 50
uLUTColorBar.ScalarBarLength = 0.3

# get color legend/bar for p_rghLUT in view renderView1
p_rghLUTColorBar = GetScalarBar(p_rghLUT, renderView1)
p_rghLUTColorBar.Title = '\xce\x94p, [pa]'
p_rghLUTColorBar.ComponentTitle = ''
p_rghLUTColorBar.TitleFontSize = 25
p_rghLUTColorBar.LabelFontSize = 25
p_rghLUTColorBar.AutomaticLabelFormat = 0
p_rghLUTColorBar.LabelFormat = '%-#6.1f'
p_rghLUTColorBar.RangeLabelFormat = '%-#6.1f'
p_rghLUTColorBar.ScalarBarThickness = 50
p_rghLUTColorBar.ScalarBarLength = 0.3

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(case1foam)
# ----------------------------------------------------------------
