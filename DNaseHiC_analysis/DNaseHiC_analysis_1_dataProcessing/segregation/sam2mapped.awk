#!/bin/awk -f

BEGIN{OFS="\t";}

{
  split($1,id,"#"); 
  print id[1],$2,$3,$4,genome,$10;
}
