
# Awk script to take a DAOphot catalog and write it to ds9 region file

BEGIN{
  print "# Region file format: DS9 version 4.1";
  print "# Filename: image.fits";
  print "global color=green font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source" 
#  print "linear";
}
{
    gsub(/,/," "); 
    if (RAD <= 0) { RAD=3.0/3600.; }
    if (WIDTH<=0) { WIDTH=2;}
    if (COLOR<=0) { COLOR="green";}
	
    if (NR>1) printf("fk5;circle(%s,%s,%s) # width=%s color=\"%s\" text={g=%3.1f}\n", $3, $4, RAD, WIDTH, COLOR, $6); 
}
