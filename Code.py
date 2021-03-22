from qgis.core import *
from osgeo import gdal, ogr, osr
import numpy
from qgis.utils import iface
import os.path
import processing
import subprocess
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry


import time

class Process:
    def __init__(self, ip_param, source_layer_filepath, source_layer_name, dst_layer_name):

        # Initialise variables
        self.base_contour = ip_param["bc"]
        self.contour_interval = ip_param["ci"]
        self.dst_layer_dir = ip_param["dst_dir"]
        self.show_index = ip_param["indx"]
        self.low_pass_filter = ip_param["lpf"]
        self.lpf_radius = ip_param["lpfr"]
        self.lpf_search_mode = ip_param["lpfsm"]
        self.z_factor = ip_param["zF"]
        self.z_factor_value = ip_param["zFval"]
        self.ignore_nodata = ip_param["ignd"]
        self.temp_dir = ip_param["td"]
        self.source_layer_filepath = source_layer_filepath
        self.source_layer_name = source_layer_name

        if len(dst_layer_name) == 0:
            dst_layer_name = source_layer_name + ' contour'
        self.dst_layer_name = dst_layer_name

        # Execute Process
        self.process()


    def process(self):
        """This Method Performs the processing part of the contour generation"""
        try:
            tic = time.time()

            # Creating Temporary Directory for temporarily saving intermediate files
            intermediate_dir = self.temp_dir + '/intermediate'
            if os.path.exists(intermediate_dir):
                intermediate_dir = modifyFolderName(intermediate_dir, 1)

            os.mkdir(intermediate_dir)

            # Initialize and Validate output file path
            dst_layer_url = self.dst_layer_dir + '\\' + self.dst_layer_name + '.shp'

            if os.path.exists(dst_layer_url):
                self.dst_layer_name = modifyLayerName(self.dst_layer_dir, self.dst_layer_name, 1)
                dst_layer_url = self.dst_layer_dir + '/' + self.dst_layer_name + '.shp'


            # Define Secondary Parameters
            # Shift base contour value to the elevation value of lowest contour of raster
            (minimum, maximum) = getMinMax(self.source_layer_filepath)
            mod_base_contour = minimum + (self.base_contour - minimum) % self.contour_interval

            # Obtain pixel size, spatial reference system
            layer = QgsRasterLayer(self.source_layer_filepath)
            pixel_size = layer.rasterUnitsPerPixelX()

            projection = layer.crs().toWkt()
            srs = osr.SpatialReference(wkt=projection)

            unprocessed_layer_url = intermediate_dir + "/unprocessed.shp"

            # ------------------------------------- Pre-Processing -----------------------------------------
            # Low Pass Filter
            if self.low_pass_filter is True:
                contour_input_filepath = sagaSimpleFilter(self.source_layer_filepath, self.lpf_search_mode,
                                                          self.lpf_radius, intermediate_dir)
            else:
                contour_input_filepath = self.source_layer_filepath

            # Apply Z-Factor
            if self.z_factor is True:
                contour_input_filepath = zFactor(self.z_factor_value,
                                                 contour_input_filepath,
                                                 intermediate_dir)

            # -------------------------------------Contour Extraction---------------------------------------
            parameters = {
                "INPUT": contour_input_filepath,
                "BAND": 1,
                "INTERVAL": self.contour_interval,
                "FIELD_NAME": 'Elev',
                "OFFSET": mod_base_contour,
                "OUTPUT": unprocessed_layer_url
            }
            if self.ignore_nodata is False:
                parameters["OUTPUT"] = dst_layer_url

            processing.run('gdal:contour', parameters)

            # ---------------------------------------- Post-Processing ----------------------------------------
            # Remove Boundary errors due to no-data region
            if self.ignore_nodata is True:
                fixBoundary(self.source_layer_filepath, intermediate_dir,
                            parameters['OUTPUT'],dst_layer_url, pixel_size,
                            srs, self.lpf_radius)

            # Load Layer to QGIS Workspace
            output_layer = iface.addVectorLayer(dst_layer_url, '', 'ogr')

            # Display Index Contour Lines
            if self.show_index is True:
                showIndex(output_layer, self.base_contour, self.contour_interval)

            toc = time.time()
            print("Time Taken: ", toc-tic)

        except:         # Error Handling
            print("An error has occurred")
            return
        return



def extractSourceLayerName(filepath):
    """
        Extract Layer name from filepath
        :Exammple:   if  filepath = "D:/folder/name.ext"
                    then layer name = "name"

        :param filepath: Path to the Layer
        :type filepath: string

        :returns layer_name: Extracted Name
        :rtype: string
    """
    filepath = filepath[::-1]   # Reverse string
    i = 0
    j = 0

    # Obtain positions of first '/' and '.'
    for c in filepath:
        i = i + 1
        if c == '.':
            j = i
        if c == '/' or c=='\\':
            i = i - 1
            break

    layer_name = filepath[j:i]      # Obtain string between '/' and '.'
    layer_name = layer_name[::-1]   # Reverse string
    return layer_name


def getMinMax(raster_layer_source):
    """
        Returns Min and Max value of a raster

        :param raster_layer_source: Path to the Layer
        :type raster_layer_source: string

        :returns minimum, maximum: Min and Max Value of raster
        :rtype: Integer Tuple
    """
    layer = QgsRasterLayer(raster_layer_source)
    extent = layer.extent()
    stats = layer.dataProvider().bandStatistics(1, QgsRasterBandStats.All, extent, 0)
    minimum = stats.minimumValue
    maximum = stats.maximumValue
    return minimum, maximum


def modifyFolderName(folder_path, i=1):
    """
        Returns a unique folder name in parent directory
        by modifying input folder name
        Example:    if folder_path = "../existingfolder"
                    then folder name is modified to "../existingfolder_xx"
                    Here xx is the first number, which gives a unique
                    name in the parent directory

        :param folder_path: Path of the folder
        :type folder_path: string

        :param i: Number added to the new name. Pass default value of 1
        :type i: int

        :returns modified_name: Path of the modified folder name
        :rtype: string
    """
    # Add number after folder name
    modified_name = folder_path + '_' + str(i)

    # If the new name is unique in parent dir.,then return the name
    if not os.path.exists(modified_name):
        return modified_name

    # If new name is not unique in parent dir.,
    # repeat the process with the next number
    return modifyFolderName(folder_path, i + 1)


def modifyLayerName(layer_dir, layer_name, i=1):
    """
        Returns a unique layer name in parent directory
        by modifying input layer name
        Example:    if layer name = "existinglayer"
                    then layer name is modified to "existinglayer_xx"
                    Here xx is the first number, which gives a unique
                    layer name in the parent directory

        Note: The extention is only limited to .shp for this function

        :param layer_dir: Path of the parent folder
        :type layer_dir: string

        :param layer_name: Name of the Layer
        :type layer_name: string

        :param i: Number added to the new name. Pass default value of 1
        :type i: int

        :returns mod_layer_name: Path of the modified layer name
        :rtype: string
    """
    # Add number after layer name
    mod_layer_name = layer_name + '_' + str(i)
    path = layer_dir + '/' + mod_layer_name + '.shp'

    # If the new name is unique in parent dir., then return the name
    if not os.path.exists(path):
        return mod_layer_name

    # If new name is not unique in parent dir.,
    # repeat the process with the next number
    return modifyLayerName(layer_dir, layer_name, i+1)


def createVectorLayer(name, dst_dir, srs, feature):
    """
        Creates a Vector Layer,

        :param name: Name of the Layer
        :type name: string

        :param dst_dir: Path of Folder for saving Layer
        :type dst_dir: string

        :param srs: Spatial Reference System for the layer
        :type srs: osgeo.osr.SpatialReference

        :param feature: Enum of Feature type
        :type feature: int
                        Available values:
                        - 1: Point feature
                        - 2: MultiLineString feature
                        - 3: Polygon feature

        :returns path: Path of the created vector layer
        :rtype: string
    """
    # Define Feature enum into feat
    feat = {
        1: ogr.wkbPoint,
        2: ogr.wkbMultiLineString,
        3: ogr.wkbPolygon,
    }

    # Specify datasource/directory
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.CreateDataSource(dst_dir)

    # Create Vector Layer
    vlayer = datasource.CreateLayer(name, srs, feat[feature])

    # Add ID and Elev Fields
    vlayer.CreateField(ogr.FieldDefn("ID", ogr.OFTInteger))
    vlayer.CreateField(ogr.FieldDefn("Elev", ogr.OFTReal))

    path = dst_dir + '/' + name + '.shp'
    return path


def sagaSimpleFilter(input_layer_source, search_mode, radius, output_dir):
    """
        Smoothens the Raster layer by applying the Simple Filter Algorithm
        from SAGA Module.

        :param input_layer_source: Path of the input Layer
        :type input_layer_source: string

        :param search_mode: Filter Shape
        :type search_mode: string
                            Available Values (Case sensitive):
                            - 'Square'
                            - 'Circle'
        :param radius: Radius of the Filter Shape (in terms of pixels)
        :type radius: int

        :param output_dir: Path of the folder for saving output file
        :type output_dir: string

        :returns output_file: Path of filtered Raster Layer
        :rtype: string
    """

    output_file = output_dir + "\SSFtemp.sdat"  # Output of SAGA algorithm only in .sdat format

    mode = {
        'Square': 0,
        'Circle': 1
    }

    parameters = {
        'INPUT': input_layer_source,
        'MODE': mode[search_mode],
        'METHOD': 0,        # 0 corresponds to Smooth
        'RADIUS': radius,
        'RESULT': output_file
        }

    processing.run('saga:simplefilter', parameters)
    return output_file


def deleteDirectory(directory_path):
    """ Deletes the directory along with the files in it """
    # Delete Contents of Directory
    for file in os.scandir(directory_path):
        # Delete Child Directory
        if file.is_dir():
            deleteDirectory(file.__fspath__())

        # Delete Files in the Directory
        else:
            try:
                os.unlink(file.__fspath__())
            # In case of error, skip deletion for that file
            except PermissionError:
                continue

    # Delete the Directory
    try:
        os.rmdir(directory_path)
    # In case of error, skip deletion of this directory
    except OSError:
        return


def fixBoundary(rlayer_url, i_dir, input_contour_url,
                output_contour_url, pixel_size, srs, lpf_radius):
    """
        Removes the error region around the boundary between
        cells with data and cells with nodata

        :param rlayer_url: Path of the input raster layer
        :type rlayer_url: string

        :param i_dir: Path to the intermediate directory
                      All files created in this functio, except the final file,
                      are stored here
        :type i_dir: string

        :param input_contour_url: Path of the input contour layer
        :type input_contour_url: string

        :param output_contour_url: Path for final processed layer
        :type output_contour_url: string

        :param pixel_size: Size of pixel of input layer
        :type pixel_size: int/float/double

        :param srs: Spatial Reference System of input layer
        :type srs: osgeo.osr.SpatialReference

        :param lpf_radius: Radius used for saga simple filter
        :type lpf_radius: int

        :returns: None
    """
    # Polygonize raster, with no data region
    polygon_w_nodata = i_dir + "/polygon_w_nodata.shp"
    cmdline = ['gdal_polygonize.py', rlayer_url, polygon_w_nodata]
    subprocess.run(cmdline, shell=True)



    #Polygonize raster, without no data region
    polygon_wo_nodata = i_dir + "/polygon_wo_nodata.shp"
    cmdline = ['gdal_polygonize.py', '-mask', rlayer_url, rlayer_url, polygon_wo_nodata]
    #-mask parameter ensure that nodata region is ignored.
    subprocess.run(cmdline, shell=True)


    #Extracting Polygon corresponding to Nodata region
    polygon_only_nodata = createVectorLayer("polygon_only_nodata", i_dir, srs, 3)
    parameters = {
        "INPUT": polygon_w_nodata,
        "OVERLAY": polygon_wo_nodata,
        "OUTPUT": polygon_only_nodata
    }
    processing.run('native:difference', parameters)


    #Creating ~70% buffer near boundary of data and nodata region
    boundary_buffer = createVectorLayer("boundary_buffer", i_dir, srs, 3)
    parameters = {
        "INPUT": polygon_only_nodata,
        "DISTANCE": pixel_size*(0.7 + 0.07  * lpf_radius) * (lpf_radius + 1),   # Empirical Formula
        "SEGMENTS": 1,
        "END_CAP_STYLE": 2, #2 = square, 1 = flat
        "JOIN_STYLE": 1, #1 = Miter
        "MITER_LIMIT": 2,
        "DISSOLVE": True,
        "OUTPUT": boundary_buffer
    }
    processing.run('native:buffer', parameters)


    #Removing Boundary Error: Unprocessed Contour - Boundary Buffer
    parameters = {
        "INPUT": input_contour_url,
        "OVERLAY": boundary_buffer,
        "OUTPUT": output_contour_url
    }
    processing.run('native:difference', parameters)

    return


def merge(layer_list, op_dir):
    """
        Merges all layers in the layer_list into a single layer

        :param layer_list: List of layer filepaths to be merged
        :type layer_list: list of string type

        :param op_dir: Path of the folder for saving merged layer
        :type op_dir: string

        :returns  path: Path to the merged layer
        :rtype: string
    """
    parameters = {
        "INPUT": layer_list,
        "PCT": False,
        "SEPARATE": False,
        "OUTPUT": op_dir + '/Merge.tif'
    }
    processing.run("gdal:merge", parameters)
    path = [parameters["OUTPUT"]]
    return path


def zFactor(factor, ip_layer_path, op_dir):
    """
        Generates a new raster with all cell values
        increased by the input factor

        :param factor: Z-factor value to be used
        :type factor: int

        :param ip_layer_path: Path of the input layer
        :type ip_layer_path: string

        :param op_dir: Path of the folder for saving the generated raster
        :type op_dir: string

        :returns op_path: Path to the new raster
        :rtype: string
    """
    layer = QgsRasterLayer(ip_layer_path)

    op_path = op_dir + '/zFactored.tif'
    op_format = "GTIFF"
    op_extent = layer.extent()
    width = layer.width()
    height = layer.height()

    # Define Raster Calculator Entries
    entries = []
    boh1 = QgsRasterCalculatorEntry()
    boh1.ref = 'boh@1'
    boh1.raster = layer     # Set input layer as boh1 entry
    boh1.bandNumber = 1
    entries.append(boh1)

    # Perform Calculation
    calc = QgsRasterCalculator(
        'boh@1 * ' + str(factor),       # Expression used in Raster Calculator
        op_path,
        op_format,
        op_extent,
        width,
        height,
        entries
    )
    calc.processCalculation()

    return op_path


def showIndex(layer, base_contour, interval):
    """
        Function to Format Contour lines by adding Index contours

        :param layer: Contour layer in which Index contours are to be shown
        :type layer: qgis._core.QgsVectorLayer

        :param base_contour: Base Contour value as a reference for selecting index contour lines
        :type base_contour: int/float/double

        :param interval: Interval of the contour
        :type interval: int/float/double

        :returns: None
    """
    showIndexLine(layer, base_contour, interval)
    showIndexLabel(layer, base_contour, interval)
    return


def showIndexLine(layer, base_contour, interval):
    """
        This function thickens every fifth line from base contour

        :param layer: Contour layer in which Index contours are to be shown
        :type layer: qgis._core.QgsVectorLayer

        :param base_contour: Base Contour value as a reference for selecting index contour lines
        :type base_contour: int/float/double

        :param interval: Interval of the contour
        :type interval: int/float/double

        :returns: None
    """
    symbol = QgsSymbol.defaultSymbol(layer.geometryType())
    renderer = QgsRuleBasedRenderer(symbol)

    # Get the root rule
    root_rule = renderer.rootRule()

    # lighten all contour lines
    def_rule = root_rule.children()[0]
    def_rule.symbol().setOpacity(0.5)

    # create clone of the default rule
    rule = root_rule.children()[0].clone()

    # set the label name and expression
    rule.setLabel("Index")
    rule.setFilterExpression(
        f"""if( ((Elev - {str(int(base_contour))})%(5*{str(int(interval))}) = 0), Elev, Null)"""
    )
    # Rule: (Elev - base) % (5 * interval) = 0                                                      |

    # Thicken and darken lines following the rule                                                   |
    width = rule.symbol().width() * 1.5
    rule.symbol().setWidth(width)
    rule.symbol().setOpacity(1)

    # append the rule to the list of rules
    root_rule.appendChild(rule)

    # apply the renderer
    layer.setRenderer(renderer)
    layer.triggerRepaint()
    return


def showIndexLabel(layer, base_contour, interval):
    """
        This function enables labels for every fifth line from base contour

        :param layer: Contour layer in which Index contours are to be shown
        :type layer: qgis._core.QgsVectorLayer

        :param base_contour: Base Contour value as a reference for selecting index contour lines
        :type base_contour: int/float/double

        :param interval: Interval of the contour
        :type interval: int/float/double

        :returns: None
    """
    settings = QgsPalLayerSettings()

    # Enable Rule-Based Labelling
    settings.isExpression = True

    # Label Placement
    settings.placement = 2              # Line Placement
    settings.placementFlags = QgsPalLayerSettings.OnLine

    # Text Formatting
    settings.fieldName = 'Elev'
    text_format = QgsTextFormat()
    text_format.setSize(10)

    buffer = QgsTextBufferSettings()
    buffer.setEnabled(True)
    buffer.setSize(1)
    text_format.setBuffer(buffer)

    settings.setFormat(text_format)

    # Create and append a new rule
    root = QgsRuleBasedLabeling.Rule(QgsPalLayerSettings())
    rule = QgsRuleBasedLabeling.Rule(settings)
    rule.setDescription("Index")
    rule.setFilterExpression(
        f"""if( ((Elev - {str(int(base_contour))})%(5*{str(int(interval))}) = 0), Elev, Null)"""
    )
    # Rule: (Elev - base) % (5 * interval) = 0

    root.appendChild(rule)

    # Apply label configuration
    rules = QgsRuleBasedLabeling(root)
    layer.setLabeling(rules)
    layer.setLabelsEnabled(True)
    layer.triggerRepaint()

    return