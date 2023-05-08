MODULE comvert_mod

IMPLICIT NONE  

include "dimensions.h"
include "paramet.h"

REAL ap(llm+1),bp(llm+1),presnivs(llm),dpres(llm),pa,preff,	&
     	nivsigs(llm),nivsig(llm+1)
! Mars Ce qui suit vient de gcm
REAL sig(llm+1),ds(llm),aps(llm),bps(llm),pseudoalt(llm)
REAL scaleheight ! atmospheric scale height (km)

END MODULE comvert_mod
