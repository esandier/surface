Claude:

A) Bugs

* Form: for disk/annulus domains, textinputs for names of variables not OK, should be as in the rectangle case.
* When specifying surface in cylindrical coordinates, thumbnails are wrong
* Save view should also save the fact that the view is orthographic/perspective, and the thumbnail be updated accordingly.

B) immersed surfaces: 

The algorithm should be as follows:

* determine self intersection curves:
  + Determine edges of the triangulation which intersect a face, this results in a set of points, each of which belongs to an edge, and to a face. 
  + Link the points in the domain, with the rule that points which belong to the same face, or to edges of the same face,  should be linked. Each segment of the curve should be remember what are the two faces of the surface it belongs to.
* Compute intersections, **in the domain** of self-intersection curves and either contour or boundary curves. When a contour or boundary curve crosses an self-intersection curve in the domain, a change of visibility of either +1 or -1 occurs on the former, this is decided using the data about the other sheet at the intersection point, are we crossing it towards the back or the front?
* intersections of contour/boundary  with self intersection curves **in the image**, which are not intersections in the domain, do not cause a visibility change (breakpoint) **of the former** but they may cause a visibility change **of the latter**, hence induce a breakpoint on the self-intersection curve.
* Connected components will include contour/boundary/self-intersection curves.
* self-intersection curves should be plotted a bit thicker than other apparent contours. 





