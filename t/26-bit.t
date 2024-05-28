#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/setbit clrbit notbit tstbit
                              bitand bitor bitxor bitnot/;

plan tests => 0 +
              5 +    # setbit
              4 +    # clrbit
              5 +    # notbit
              4 +    # tstbit
              5 +    # bitand
              5 +    # bitor
              5 +    # bitxor
              2 +    # bitnot
              2 +    # ... extra
            + 0;

##### setbit
is_deeply([map { setbit(0,$_) } 0..16], [map { 1<<$_ } 0..16], "setbit(0,b) = 2^b");
is_deeply([map { setbit(1,$_) } 0..16], [map { (1<<$_) | 1 } 0..16], "setbit(1,b) = 2^b | 1");
is_deeply([map { setbit(13,$_) } 0..16], [map { (1<<$_) | 13 } 0..16], "setbit(13,b) = 2^b | 13");
is_deeply([map { setbit(-1,$_) } 0..16], [map { -1 } 0..16], "setbit(-1,b) = -1");
is(setbit(4096,0),4097,"setbit(4096,0) = 4097");

##### clrbit
is_deeply([map { clrbit(0,$_) } 0..16], [map { 0 } 0..16], "clrbit(0,b) = 0");
is_deeply([map { clrbit(2,$_) } 0..5], [2,0,2,2,2,2], "clrbit(2,b) = {0 for b=1, 2 otherwise}");
is(clrbit(4097,0),4096,"clrbit(4097,0) = 4096");
is(clrbit(4097,12),1,"clrbit(4097,1) = 1");

##### notbit   (combit)
is_deeply([map { notbit($_,0) } -5..5],[-6,-3,-4,-1,-2,1,0,3,2,5,4],"notbit(-5..5,0) clears lower bit");
is(notbit(4097,0),4096,"notbit(4097,0) = 4096");
is(notbit(4096,0),4097,"notbit(4096,0) = 4097");
is(notbit( 4097,15), 36865,"notbit( 4097,15) =  36865");
is(notbit(-4097,15),-36865,"notbit(-4097,15) = -36865");

##### tstbit
is(tstbit(7,0),1,"tstbit(7,0) = 1");
is(tstbit(7,2),1,"tstbit(7,2) = 1");
is(tstbit(7,3),0,"tstbit(7,3) = 0");
is(tstbit(-1,10),1,"tstbit(-1,10) = 1");


##### bitand
is(bitand( 5, 3), 1,"bitand( 5, 3) =  1");
is(bitand(-5, 3), 3,"bitand(-5, 3) =  3");
is(bitand(-5,-3),-7,"bitand(-5,-3) = -7");
is(bitand(3231,333437),1053,"bitand(3231,333437) = 1053");
is(bitand("340282366920938463481821351505477763073","-340282366920938463481821351505477763073"),1,"bitand(340282366920938463481821351505477763073, -340282366920938463481821351505477763073) = 1");

##### bitor
is(bitor( 5, 3), 7,"bitor( 5, 3) =  7");
is(bitor(-5, 3),-5,"bitor(-5, 3) = -5");
is(bitor(-5,-3),-1,"bitor(-5,-3) = -1");
is(bitor(3231,333437),335615,"bitor(3231,333437) = 335615");
is(bitor("340282366920938463481821351505477763073","-340282366920938463481821351505477763073"),-1,"bitor(340282366920938463481821351505477763073, -340282366920938463481821351505477763073) = -1");

##### bitxor
is(bitxor( 5, 3), 6,"bitxor( 5, 3) =  6");
is(bitxor(-5, 3),-8,"bitxor(-5, 3) = -8");
is(bitxor(-5,-3), 6,"bitxor(-5,-3) =  6");
is(bitxor(3231,333437),334562,"bitxor(3231,333437) = 334562");
is(bitxor("340282366920938463481821351505477763073","-340282366920938463481821351505477763073"),-2,"bitxor(340282366920938463481821351505477763073, -340282366920938463481821351505477763073) = -2");

##### bitnot
is_deeply([map { bitnot($_) } -9 .. 9], [8,7,6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10], "bitnot(-9 .. 9)");
is(bitnot("36893488147419103231"), "-36893488147419103232","bitnot(36893488147419103231) = -36893488147419103232");

##### extra
is_deeply([map { setbit($_,7) } -16..16], [map { bitor($_,2**7) } -16..16], "setbit(-16..16,7) = bitor(-16..16,2^7)");
is_deeply([map { clrbit($_,7) } -16..16], [map { bitand($_,bitnot(2**7)) } -16..16], "clrbit(-16..16,7) = bitand(-16..16,bitnot(2^7))");
