#!/bin/csh 
# -i in place BUT save bu file with extension .sed
#generically 's/old/new/'    \ = delimitting char
set old = _r_2 
set new = 
set FILE = *90
#so need to adjust search and replace strings here
#sed -i .sed 's/include \"include/include \"..\/..\/include/' cable_*
#sed -i .sed 's/old              /new                      /' cable_*

echo $old $new $FILE

sed -i .sed 's/'$old'/'$new'/'  $FILE

