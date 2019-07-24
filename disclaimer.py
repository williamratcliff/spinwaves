"""
Disclaimer
==========

This software was developed at the National Institute of Standards and Technology at the NIST Center for Neutron Research by employees of the Federal Government in the course of their official duties. Pursuant to title 17 section 105* of the United States Code this software is not subject to copyright protection and is in the public domain. The SPINAL software package is an experimental spinwave analysis system. NIST assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. The use of certain trade names or commercial products does not imply any endorsement of a particular product, nor does it imply that the named product is necessarily the best product for the stated purpose. We would appreciate acknowledgment if the software is used.

*Subject matter of copyright: United States Government works

Copyright protection under this title is not available for any work of the United States Government, but the United States Government is not precluded from receiving and holding copyrights transferred to it by assignment, bequest, or otherwise."""


"""This script will add a disclaimer to all python files in the same directory
or any subdirectories."""

import os
import string

disclaimer = """\"\"\"
Disclaimer
==========

This software was developed at the National Institute of Standards and Technology at the NIST Center for Neutron Research by employees of the Federal Government in the course of their official duties. Pursuant to title 17 section 105* of the United States Code this software is not subject to copyright protection and is in the public domain. The SPINAL software package is an experimental spinwave analysis system. NIST assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. The use of certain trade names or commercial products does not imply any endorsement of a particular product, nor does it imply that the named product is necessarily the best product for the stated purpose. We would appreciate acknowledgment if the software is used.

*Subject matter of copyright: United States Government works

Copyright protection under this title is not available for any work of the United States Government, but the United States Government is not precluded from receiving and holding copyrights transferred to it by assignment, bequest, or otherwise.\"\"\"


"""

def add_disclaimer(file):
    """file is the absolute path of the python file which is having the disclaimer added.  This
    function will add disclaimer if the file's contents does not already begin with the contents
    of disclaimer.  If the disclaimer text is changed, than this functionality must be changed
    too(Perhaps just compare the disclaimer title)."""
    handle = open(file, "r+")
    file_contents = handle.read()
    if file_contents[0:len(disclaimer)] != disclaimer:
        handle.seek(0)
        handle.write(disclaimer)
        handle.write(file_contents)
        handle.close()
                
def handle_dir(dir):
    contents = os.listdir(dir)
    for item in contents:
        item = os.path.join(dir, item)
        if os.path.isdir(item):
            handle_dir(item)
        else:
            if os.path.isfile(item):
                file = os.path.split(item)[1]
                parts = string.split(file, ".")
                if len(parts) > 1 and parts[len(parts)-1] == "py":#It's a python file.
                    print item + " handled."
                    add_disclaimer(item)

                              
root_dir = os.path.dirname(__file__)
handle_dir(root_dir)