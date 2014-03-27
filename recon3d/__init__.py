"""
Recon3D -- Python-based routine for reconstructing   
         3D displacement (velocity) fields from multiple
         InSAR data sets acquired from varios vantage
         points 

Copyright (C) 2013   Brent M. Minchew
--------------------------------------------------------------------
GNU Licensed

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------

"""

import version
import _reconparser
import _recontools
import _reconsolver
import _reconutils

__version__ = version.version
__all__ = ['_reconparser','_recontools','_reconsolver','_reconutils']
