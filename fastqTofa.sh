awk '{if(NR%4==2 || NR%4==1) print $0}' ${1} | sed 's/@/>/;N'
