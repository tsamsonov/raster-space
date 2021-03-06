# -*- coding: utf-8 -*-

"""
/***************************************************************************
 RasterSpace
                                 A QGIS plugin
 Estimation of the local space size using raster algebra
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2019-07-26
        copyright            : (C) 2020 by Timofey Samsonov, Lomonosov MSU Faculty of Geography
        email                : tsamsonov@geogr.msu.ru
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Timofey Samsonov, Lomonosov MSU Faculty of Geography'
__date__ = '2020-07-26'
__copyright__ = '(C) 2020 by Timofey Samsonov, Lomonosov MSU Faculty of Geography'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

import os
import sys
import inspect

from qgis.core import QgsApplication
from .raster_space_provider import RasterSpaceProvider

class RasterSpacePlugin(object):

    def __init__(self):
        self.provider = None

    def initProcessing(self):
        """Init Processing provider for QGIS >= 3.8."""
        self.provider = RasterSpaceProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
