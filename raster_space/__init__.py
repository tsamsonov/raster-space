# -*- coding: utf-8 -*-
"""
/***************************************************************************
 RasterSpace
                                 A QGIS plugin
 Estimation of the local space size using raster algebra
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2019-07-26
        copyright            : (C) 2019 by Timofey Samsonov, Lomonosov MSU Faculty of Geography
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
 This script initializes the plugin, making it known to QGIS.
"""

__author__ = 'Timofey Samsonov, Lomonosov MSU Faculty of Geography'
__date__ = '2019-07-26'
__copyright__ = '(C) 2019 by Timofey Samsonov, Lomonosov MSU Faculty of Geography'


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load RasterSpace class from file RasterSpace.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .raster_space import RasterSpacePlugin
    return RasterSpacePlugin()