import sys
sys.path.append('/home/bst/other/ehalper/python/pysam-0.5/build/lib.linux-x86_64-2.6')
import pysam
samuse=sys.argv[1]
samwrite=sys.argv[2]
samfile=pysam.Samfile(samuse,'rb')
reads=pysam.Samfile(samwrite, "wb", template=samfile)

def get_bit(byteval,idx):
    return ((byteval&(1<<idx))!=0);

for alignedread in samfile.fetch():
    reads.write(alignedread)
    originalflag=alignedread.flag
    if(alignedread.tags[0][0] == 'X0'):
        try:   
            for newrec in alignedread.opt('XA').split(';'):
                if(newrec != ''):
                    s=newrec.split(',')
                    a = pysam.AlignedRead()
                    a = alignedread
                    a.flag=originalflag
                    if(get_bit(alignedread.flag,4) and int(s[1]) < 0 and get_bit(alignedread.flag,7)):
                        a.flag = int(0x80)
                    elif(get_bit(alignedread.flag,4) and int(s[1]) > 0 and get_bit(alignedread.flag,7)):
                        a.flag= int(0x10) + int(0x80)
                    elif(get_bit(alignedread.flag,4) and int(s[1]) < 0 and get_bit(alignedread.flag,6)):
                        a.flag= int(0x40)
                    elif(get_bit(alignedread.flag,4) and int(s[1]) > 0 and get_bit(alignedread.flag,6)):
                        a.flag=int(0x10) + int(0x40)
                    elif(int(s[1]) > 0 and get_bit(alignedread.flag,7)):
                        a.flag= int(0x80)     
                    elif(int(s[1]) < 0 and get_bit(alignedread.flag,7)):
                        a.flag= int(0x10) + int(0x80)
                    elif(int(s[1]) > 0 and get_bit(alignedread.flag,6)):
                        a.flag= int(0x40)
                    elif(int(s[1]) < 0 and get_bit(alignedread.flag,6)):
                        a.flag=int(0x10) + int(0x40)
                    a.pos=abs(int(s[1]))-1
                    a.flag+=1 ## Set first bit to READ IS PAIRED
                    if(s[0] == 'Y'):
                        rname=(24-1)
                    elif(s[0]== 'X'):
                        rname=(23-1)
                    else:
                        rname=int(s[0])-1  ## WE NEED THE MINUS 1 BECAUSE PYSAM TRYS TO CORRECT ZERO OFFSET WHEN WRITING

                    a.rname=rname

                    #print(s[3])
                    a.tags=[]
                    a.tags = [('NM',int(s[3]))]
                    reads.write(a)
        except:
            pass
                    
        

reads.close()
samfile.close()
