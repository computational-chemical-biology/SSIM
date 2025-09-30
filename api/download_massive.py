import ftplib
import re
import sys
import os

#ftp://massive-ftp.ucsd.edu/v02/MSV000083183/
def downloadMassIVE(massive_url):
    ftp = ftplib.FTP("massive-ftp.ucsd.edu")
    ftp.login('anonymous', 'password')
    dr = re.sub('.+edu/(v\d+/M.+)', '\\1', massive_url)
    ftp.cwd(dr)
    # print mzxml files 
    #fls = ftp.nlst('ccms_peak/')
    fls = ftp.nlst('raw/')
    if len(fls):
        for f in fls:
            fo = f.split('/')[-1]
            print(f'Downloading {fo}...')
            fo = fo.replace('mzXML', 'mzML')
            os.system(f'wget -O {fo} https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000085496%2Fccms_peak%2F{fo}')
            #ftp.retrbinary("RETR " +f, open(fo, 'wb').write)
    else:
        print('No mzxml file found')
    ftp.quit()


if __name__=='__main__':
    downloadMassIVE(sys.argv[1])
