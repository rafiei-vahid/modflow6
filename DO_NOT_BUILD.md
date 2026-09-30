# DO NOT BUILD

This branch is a DELIBERATELY BROKEN control build: Q1041 E v3 control 2 (matrix).
Its parent commit 1de802670fcd443b185cd6d9d2d97398350e8ea7 breaks the matrix half of DRY_CELL_ROUTE_MASS on purpose,
so that the acceptance test for the real patch (branch q1041-e-v3,
e3f1098cb634d7b68c6e32976b521c3b0d198999) could prove it detects that failure.

Never build it for a real run, never deploy it, never cite a result from it.
The real change is on q1041-e-v3.
