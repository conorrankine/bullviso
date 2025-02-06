"""
BULLVISO
Copyright (C) 2024  Conor D. Rankine

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either Version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

###############################################################################
################################### CLASSES ###################################
###############################################################################

class BVBarcode:

    def __init__(self, barcode: tuple, grouped_barcode: tuple = None):
        """
        Instantiates a bullvalene isomer barcode.

        Args:
            barcode (tuple): A tuple of 10 integers in the range 0-9 inclusive
                representing the bullvalene isomer barcode.

        Raises:
            ValueError: If the barcode i) is not of length 10, ii) contains
                non-integer elements, or iii) contains integers less than 0
                or greater than 9 (i.e. not in the range 0-9 inclusive).
        """

        if len(barcode) == 10 and all(
            isinstance(x, int) and 0 <= x <= 9 for x in barcode
        ):        
            self._barcode = barcode
        else:
            raise ValueError(
                f'{barcode} is not a valid barcode'
            )
        
        if grouped_barcode is not None:
            if len(grouped_barcode) == 10 and all(
                isinstance(x, int) and 0 <= x <= 9 for x in grouped_barcode
            ):        
                self._grouped_barcode = grouped_barcode
            else:
                raise ValueError(
                    f'{grouped_barcode} is not a valid grouped barcode'
                )
        else:
            self._grouped_barcode = barcode

    @property
    def barcode(self) -> tuple:
        """
        Returns the bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The bullvalene isomer barcode as a tuple.
        """

        return self._barcode
    
    @property
    def grouped_barcode(self) -> tuple:
        """
        Returns the grouped bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The grouped bullvalene isomer barcode as a tuple.
        """

        return self._grouped_barcode
    