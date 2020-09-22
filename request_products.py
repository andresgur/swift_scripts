# @Author: Andrés Gúrpide <agurpide>
# @Date:   17-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 17-09-2020

from swifttools.xrt_prods import XRTProductRequest

myReq = XRTProductRequest('agurpidelash@irap.omp.eu', silent=False)
myReq.setGlobalPars(centroid=True, centMeth="simple", posErr=1, name="Holmberg II X-1")

myReq.addLightCurve(binMeth='counts', pcCounts=20, wtCounts=30, dynamic=True)
