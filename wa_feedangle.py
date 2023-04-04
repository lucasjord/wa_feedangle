
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WAIPSUVData

import numpy as np, argparse, sys, os, datetime
from numpy import cos, sin, tan, arcsin, arccos, arctan2

from astropy.coordinates import SkyCoord
from astropy.time import Time as aTime
import astropy.units as u

################################################

def main():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument('inname',
                        help='Catalogue INNAME',type=str)
    parser.add_argument('aips_id',
                        help='AIPS ID, e.g 666',type=int)
    # optional arguments
    parser.add_argument('-k','--klass',
                        help='Catalogue KLASS. Default UVDATA',type=str,default='UVDATA')
    parser.add_argument('-q','--seq',
                        help='Catalogue SEQ. Default 1',type=int,default=1)
    parser.add_argument('-d','--disk',
                        help='Catalogue DISK. Default 1',type=int,default=1)
    parser.add_argument('-i','--cl_in',
                        help='CL table in. Default HIGHVER',type=int,default=None)    
    parser.add_argument('-o','--cl_out',
                        help='CL table out. Default HIGHVER+1',type=int,default=None)
    parser.add_argument('-s','--sn_out',
                        help='SN table out. Default HIGHVER+1',type=int,default=None)
    parser.add_argument('--pang',
                        help='Also add parallactic angle correction.',
                        action="store_true",default=False)
    parser.add_argument('--swpol',
                        help='Swap correction sign.',
                        action="store_true",default=False)
    args = parser.parse_args()
    '''
    Look for input catalogue.
    '''
    AIPS.userno = args.aips_id

    indata = AIPSUVData(args.inname,args.klass,args.disk,args.seq)
    if not indata.exists():
        sys.exit('{}.{}.{}.{} does not exist on AIPSID {}'.format(args.inname,
                                    args.klass,args.disk,args.seq,args.aips_id))
    '''
    Check if CL tables are okay.
    '''
    cl_tables = [s[0] for s in indata.tables if 'CL' in s[-1]]
    #
    if args.cl_in==None and args.cl_out==None:
        args.cl_in  = max(cl_tables)
        args.cl_out = args.cl_in + 1
    elif args.cl_in!=None and args.cl_out==None:
        args.cl_out = args.cl_in + 1
    #
    cl_high = indata.table_highver
    if args.cl_in > cl_high:
        sys.exit('CL{} does not exist. Cannot use to calibrate data. Highest version is CL{}'.format(args.cl_in,cl_high))
    '''
    Running feed angle correction.
    '''
    sn_tables = [s[0] for s in indata.tables if 'SN' in s[-1]] + [0]
    sn_high   = max(sn_tables)
    if args.sn_out==None: 
        args.sn_out = sn_high + 1
    elif sn_high>=args.sn_out:
        print 'Deleting old SN tables > {}'.format(args.sn_out)
        while indata.table_highver('SN')>=args.sn_out:
            indata.zap_table('SN',indata.table_highver('SN'))
    #
    success   = run_wa_pang(indata,1,args.swpol,args.cl_in,args.cl_out,args.sn_out)
    if success==0:
        print 'WAPANG did not run.' 
    if success==1:
        print 'Task appears to have ended successfully.'

################################################

def runtbin(outdata,intext):
    tbin          = AIPSTask('TBIN')
    tbin.outdata  = outdata
    tbin.intext   = intext
    tbin()

def runclcal(indata,snver,cl_in,cl_out,refant,interpol='',dobtween=1):
    clcal          = AIPSTask('CLCAL')
    clcal.indata   = indata
    clcal.refant   = refant
    clcal.snver    = snver
    clcal.inver    = 0
    clcal.gainver  = cl_in
    clcal.gainuse  = cl_out
    clcal.interpol = interpol
    clcal.dobtween = dobtween
    clcal()

class telescope:
    def __init__(self,name,lat,long):
        self.lat  = lat   # deg
        self.long = long  # deg
        self.name = name

def run_wa_pang(indata,refant,swpol,pang,cl_in=3,cl_out=4,sn_out=0):
    '''
    Apply feed rotation corrections for Warkworth30m being 
    beam-waveguide/Naysmith + 6x'parascope' mirrors 
    - LJH 2021/11/30 
    '''
    if len(indata.stokes)==1:
        print '###################################'
        print 'ANT only has a single pol, skipping'
        print '###################################'
        return 0
    wa_num  = indata.antennas.index('WA') + 1
    wa      = telescope('WA',-36.43316,174.66295) # currently hardcoded
    L       = wa.lat*u.deg.to(u.rad)  # this step incase needs to be changed to another antenna
    wizdata = WAIPSUVData(indata.name,indata.klass,indata.disk,indata.seq)
    
    ## make SU tables and get visibilities with Wizardry
    print 'WAPANG: Getting source RA/DECs.'
    su = [{s.id__no:[s.raapp,s.decapp]} for s in indata.table('SU',1)]
    SU = {}
    for s in su: SU.update(s) #making SU dictionary
    
    print 'WAPANG: Getting visibilities using wizardry.'
    vis = []
    for v in wizdata:
        if v.baseline[0]==v.baseline[1]:
            vis += [[v.baseline[0],v.time]+SU[v.source]+[v.source,v.inttim]]
    V = np.array(vis)
    t, ra, dec, sid, dt = V[V[:,0]==wa_num,1:].T.tolist() 
    # time (fday), ra (deg), dec (deg), sourceID, int time (sec)
    
    ## getting start time and convert to LST
    start = aTime(datetime.datetime.strptime(indata.header['date_obs'],'%Y-%m-%d'))
    T = start + t*u.day
    T.delta_ut1_utc = 0. # doing this to suppress errors. Might be bad #yolo
    print 'WAPANG: Calculating WA sidereal time.'
    st = T.sidereal_time('mean',longitude=wa.long*u.deg).value
    
    ## calculate hour angle
    print 'WAPANG: Calculating hour angle.'
    ha_rad = (st*u.hourangle).to(u.rad)-(ra*u.deg).to(u.rad)
    
    ## calculate elevation angle
    print 'WAPANG: Calculating elevation angle.'
    sine = sin(L)*sin((dec*u.deg).to(u.rad)) + cos(L)*cos((dec*u.deg).to(u.rad))*cos(ha_rad)
    el   = arcsin(sine)

    ## calculate azimuth angle
    print 'WAPANG: Calculating azimuth angle.'
    tan_an = -cos((dec*u.deg).to(u.rad))*sin(ha_rad) 
    tan_ad = ( cos(L)*sin((dec*u.deg).to(u.rad))
             - sin(L)*cos((dec*u.deg).to(u.rad))*cos(ha_rad) )
    az     = arctan2(tan_an,tan_ad)

    ## calculate parallactic angle 
    print 'WAPANG: Calculating parrallactic angle.'
    tan_pn = cos(L)*sin(ha_rad)
    tan_pd = ( sin(L)*cos((dec*u.deg).to(u.rad)) 
             - cos(L)*sin((dec*u.deg).to(u.rad))*cos(ha_rad) )
    pa = arctan2(tan_pn,tan_pd)

    ## creating correction factor assuming that PANG has already been run
    ## pang adds/subtracts parallactic angle to antenna rcp/lcp
    print 'WAPANG: Calculating feed angle.'
    if pang==False: feed_angle = (       az - el ).to(u.rad).value*(-1)**swpol
    else:           feed_angle = ( -pa + az - el ).to(u.deg).value*(-1)**swpol  # if no PANG done first

    ## going to make SN table to add to Warkworth phase
    header  = ["XTENSION= 'BINTABLE'           / extension type",
       'BITPIX  =                    8 / printable ASCII codes',
       'NAXIS   =                    2 / Table is a matrix',
       'NAXIS1  =                  999 / Max. no. of characters/pass',
       'NAXIS2  =                  732 / Number of entries in table',
       'PCOUNT  =                    0 / Random parameter count',
       'GCOUNT  =                    1 / Group count',
       'NOPASS  =                    1 / Number of passes thru table',
       'TFIELDS =                   26 / Number of fields in each row',
       "EXTNAME = 'AIPS SN '           / AIPS table file",
       'EXTVER  =                    1 / Version Number of table',
       'TBCOL1  =                    9 / Starting char. pos. of field',
       "TFORM1  = 'D24.15  '           / Fortran format of field  1",
       'TFDIM1  =                    1 / Dimension of field  1',
       "TTYPE1  = 'TIME                    '           / type (heading) of field  1",
       "TUNIT1  = 'DAYS    '           / physical units of field  1",
       'TBCOL2  =                   33 / Starting char. pos. of field',
       "TFORM2  = 'E15.6   '           / Fortran format of field  2",
       'TFDIM2  =                    1 / Dimension of field  2',
       "TTYPE2  = 'TIME INTERVAL           '           / type (heading) of field  2",
       "TUNIT2  = 'DAYS    '           / physical units of field  2",
       'TBCOL3  =                   48 / Starting char. pos. of field',
       "TFORM3  = 'I11     '           / Fortran format of field  3",
       'TFDIM3  =                    1 / Dimension of field  3',
       "TTYPE3  = 'SOURCE ID               '           / type (heading) of field  3",
       "TUNIT3  = '        '           / physical units of field  3",
       'TBCOL4  =                   59 / Starting char. pos. of field',
       "TFORM4  = 'I11     '           / Fortran format of field  4",
       'TFDIM4  =                    1 / Dimension of field  4',
       "TTYPE4  = 'ANTENNA NO.             '           / type (heading) of field  4",
       "TUNIT4  = '        '           / physical units of field  4",
       'TBCOL5  =                   70 / Starting char. pos. of field',
       "TFORM5  = 'I11     '           / Fortran format of field  5",
       'TFDIM5  =                    1 / Dimension of field  5',
       "TTYPE5  = 'SUBARRAY                '           / type (heading) of field  5",
       "TUNIT5  = '        '           / physical units of field  5",
       'TBCOL6  =                   81 / Starting char. pos. of field',
       "TFORM6  = 'I11     '           / Fortran format of field  6",
       'TFDIM6  =                    1 / Dimension of field  6',
       "TTYPE6  = 'FREQ ID                 '           / type (heading) of field  6",
       "TUNIT6  = '        '           / physical units of field  6",
       'TBCOL7  =                   92 / Starting char. pos. of field',
       "TFORM7  = 'E15.6   '           / Fortran format of field  7",
       'TFDIM7  =                    1 / Dimension of field  7',
       "TTYPE7  = 'I.FAR.ROT               '           / type (heading) of field  7",
       "TUNIT7  = 'RAD/M**2'           / physical units of field  7",
       'TBCOL8  =                  107 / Starting char. pos. of field',
       "TFORM8  = 'I11     '           / Fortran format of field  8",
       'TFDIM8  =                    1 / Dimension of field  8',
       "TTYPE8  = 'NODE NO.                '           / type (heading) of field  8",
       "TUNIT8  = '        '           / physical units of field  8",
       'TBCOL9  =                  118 / Starting char. pos. of field',
       "TFORM9  = 'E15.6   '           / Fortran format of field  9",
       'TFDIM9  =                    1 / Dimension of field  9',
       "TTYPE9  = 'MBDELAY1                '           / type (heading) of field  9",
       "TUNIT9  = 'SECONDS '           / physical units of field  9",
       'TBCOL10 =                  133 / Starting char. pos. of field',
       "TFORM10 = 'E15.6   '           / Fortran format of field 10",
       'TFDIM10 =                    1 / Dimension of field 10',
       "TTYPE10 = 'DISP 1                  '           / type (heading) of field 10",
       "TUNIT10 = 'SEC/M**2'           / physical units of field 10",
       'TBCOL11 =                  148 / Starting char. pos. of field',
       "TFORM11 = 'E15.6   '           / Fortran format of field 11",
       'TFDIM11 =                    1 / Dimension of field 11',
       "TTYPE11 = 'DDISP 1                 '           / type (heading) of field 11",
       "TUNIT11 = 'S/S/M**2'           / physical units of field 11",
       'TBCOL12 =                  163 / Starting char. pos. of field',
       "TFORM12 = 'E15.6   '           / Fortran format of field 12",
       'TFDIM12 =                    8 / Dimension of field 12',
       "TTYPE12 = 'REAL1                   '           / type (heading) of field 12",
       "TUNIT12 = '        '           / physical units of field 12",
       'TBCOL13 =                  178 / Starting char. pos. of field',
       "TFORM13 = 'E15.6   '           / Fortran format of field 13",
       'TFDIM13 =                    8 / Dimension of field 13',
       "TTYPE13 = 'IMAG1                   '           / type (heading) of field 13",
       "TUNIT13 = '        '           / physical units of field 13",
       'TBCOL14 =                  193 / Starting char. pos. of field',
       "TFORM14 = 'E15.6   '           / Fortran format of field 14",
       'TFDIM14 =                    8 / Dimension of field 14',
       "TTYPE14 = 'DELAY 1                 '           / type (heading) of field 14",
       "TUNIT14 = 'SECONDS '           / physical units of field 14",
       'TBCOL15 =                  208 / Starting char. pos. of field',
       "TFORM15 = 'E15.6   '           / Fortran format of field 15",
       'TFDIM15 =                    8 / Dimension of field 15',
       "TTYPE15 = 'RATE 1                  '           / type (heading) of field 15",
       "TUNIT15 = 'SEC/SEC '           / physical units of field 15",
       'TBCOL16 =                  223 / Starting char. pos. of field',
       "TFORM16 = 'E15.6   '           / Fortran format of field 16",
       'TFDIM16 =                    8 / Dimension of field 16',
       "TTYPE16 = 'WEIGHT 1                '           / type (heading) of field 16",
       "TUNIT16 = '        '           / physical units of field 16",
       'TBCOL17 =                  238 / Starting char. pos. of field',
       "TFORM17 = 'I11     '           / Fortran format of field 17",
       'TFDIM17 =                    8 / Dimension of field 17',
       "TTYPE17 = 'REFANT 1                '           / type (heading) of field 17",
       "TUNIT17 = '        '           / physical units of field 17",
       'TBCOL18 =                  249 / Starting char. pos. of field',
       "TFORM18 = 'E15.6   '           / Fortran format of field 18",
       'TFDIM18 =                    1 / Dimension of field 18',
       "TTYPE18 = 'MBDELAY2                '           / type (heading) of field 18",
       "TUNIT18 = 'SECONDS '           / physical units of field 18",
       'TBCOL19 =                  264 / Starting char. pos. of field',
       "TFORM19 = 'E15.6   '           / Fortran format of field 19",
       'TFDIM19 =                    1 / Dimension of field 19',
       "TTYPE19 = 'DISP 2                  '           / type (heading) of field 19",
       "TUNIT19 = 'SEC/M**2'           / physical units of field 19",
       'TBCOL20 =                  279 / Starting char. pos. of field',
       "TFORM20 = 'E15.6   '           / Fortran format of field 20",
       'TFDIM20 =                    1 / Dimension of field 20',
       "TTYPE20 = 'DDISP 2                 '           / type (heading) of field 20",
       "TUNIT20 = 'S/S/M**2'           / physical units of field 20",
       'TBCOL21 =                  294 / Starting char. pos. of field',
       "TFORM21 = 'E15.6   '           / Fortran format of field 21",
       'TFDIM21 =                    8 / Dimension of field 21',
       "TTYPE21 = 'REAL2                   '           / type (heading) of field 21",
       "TUNIT21 = '        '           / physical units of field 21",
       'TBCOL22 =                  309 / Starting char. pos. of field',
       "TFORM22 = 'E15.6   '           / Fortran format of field 22",
       'TFDIM22 =                    8 / Dimension of field 22',
       "TTYPE22 = 'IMAG2                   '           / type (heading) of field 22",
       "TUNIT22 = '        '           / physical units of field 22",
       'TBCOL23 =                  324 / Starting char. pos. of field',
       "TFORM23 = 'E15.6   '           / Fortran format of field 23",
       'TFDIM23 =                    8 / Dimension of field 23',
       "TTYPE23 = 'DELAY 2                 '           / type (heading) of field 23",
       "TUNIT23 = 'SECONDS '           / physical units of field 23",
       'TBCOL24 =                  339 / Starting char. pos. of field',
       "TFORM24 = 'E15.6   '           / Fortran format of field 24",
       'TFDIM24 =                    8 / Dimension of field 24',
       "TTYPE24 = 'RATE 2                  '           / type (heading) of field 24",
       "TUNIT24 = 'SEC/SEC '           / physical units of field 24",
       'TBCOL25 =                  354 / Starting char. pos. of field',
       "TFORM25 = 'E15.6   '           / Fortran format of field 25",
       'TFDIM25 =                    8 / Dimension of field 25',
       "TTYPE25 = 'WEIGHT 2                '           / type (heading) of field 25",
       "TUNIT25 = '        '           / physical units of field 25",
       'TBCOL26 =                  369 / Starting char. pos. of field',
       "TFORM26 = 'I11     '           / Fortran format of field 26",
       'TFDIM26 =                    8 / Dimension of field 26',
       "TTYPE26 = 'REFANT 2                '           / type (heading) of field 26",
       "TUNIT26 = '        '           / physical units of field 26",
       'NO_ANT  =            4', 'NO_POL  =            2',
       'NO_IF   =            8', 'NO_NODES=            0',
       'MGMOD   =   1.00000000000000000D+00',
       'APPLIED =                    F', 'REVISION=            1',
       'SNORIGIN=            0',
       "HISTORY TBOUT  /SN table version  1 of INNAME='S001K_C     .UVDATA.   1'",
       'HISTORY TBOUT  / copied to the text file HOME:sn_poo.TBOUT',
       'HISTORY TBOUT  / 09-DEC-2021 13:19:05', 'END',
       'COL. NO.      1                       2              3          4          5          6          7              8          9             10             11             12             13             14             15             16             17         18             19             20             21             22             23             24             25             26',
       'ROW   TIME                    TIME INT       SOURCE I   ANTENNA    SUBARRAY   FREQ ID    I.FAR.RO       NODE NO.   MBDELAY1       DISP 1         DDISP 1        REAL1          IMAG1          DELAY 1        RATE 1         WEIGHT 1       REFANT 1   MBDELAY2       DISP 2         DDISP 2        REAL2          IMAG2          DELAY 2        RATE 2         WEIGHT 2       REFANT 2',
       'NUMBER   DAYS                    DAYS                                                       RAD/M**2                  SECONDS        SEC/M**2       S/S/M**2                                     SECONDS        SEC/SEC                                  SECONDS        SEC/M**2       S/S/M**2                                     SECONDS        SEC/SEC',
       '***BEGIN*PASS***']

    ncpx, nstk, nfrq, nif, nra, ndec, bgs = indata.header['naxis']
    
    if nif==1: header[0] = "XTENSION= 'TABLE   '           / extension type"
    else:      header[0] = "XTENSION= 'BINTABLE'           / extension type"

    header[4]   = 'NAXIS2  =              {:>7s} / Number of entries in table'.format(str(len(t)))

    header[10]  = 'EXTVER  =                   {:>2.0f} / Version Number of table'.format(sn_out)

    header[68]  = 'TFDIM12 =                   {:2.0f} / Dimension of field 12'.format(nif)
    header[73]  = 'TFDIM13 =                   {:2.0f} / Dimension of field 13'.format(nif)
    header[78]  = 'TFDIM14 =                   {:2.0f} / Dimension of field 14'.format(nif)
    header[83]  = 'TFDIM15 =                   {:2.0f} / Dimension of field 15'.format(nif)
    header[88]  = 'TFDIM16 =                   {:2.0f} / Dimension of field 16'.format(nif)
    header[93]  = 'TFDIM17 =                   {:2.0f} / Dimension of field 17'.format(nif)
    header[113] = 'TFDIM21 =                   {:2.0f} / Dimension of field 21'.format(nif)
    header[118] = 'TFDIM22 =                   {:2.0f} / Dimension of field 22'.format(nif)
    header[123] = 'TFDIM23 =                   {:2.0f} / Dimension of field 23'.format(nif)
    header[128] = 'TFDIM24 =                   {:2.0f} / Dimension of field 24'.format(nif)
    header[133] = 'TFDIM25 =                   {:2.0f} / Dimension of field 25'.format(nif)
    header[138] = 'TFDIM26 =                   {:2.0f} / Dimension of field 26'.format(nif)

    header[141] = 'NO_ANT  =           {:>2.0f}'.format(len(indata.antennas))
    header[143] = 'NO_IF   =           {:>2.0f}'.format(nif)

    ## print solutions to file to make SN_OUT
    if os.path.exists('wa_pang.TBIN'): os.remove('wa_pang.TBIN')
    with open('wa_pang.TBIN','w+') as w:
        for row in header: print >> w, row
        for j in range(len(t)):
            for i in range(nif):
                if i==0: 
                    print >> w, '{0:8.0f}{1:>24.15E}{2:>15.6E}{3:>11.0f}{4:>11.0f}{5:>11.0f}{6:>11.0f}{7:>15.3f}{8:>11.0f}{9:>15.3f}{10:>15.3f}{11:>15.3f}{12:>15.6f}{13:>15.6f}{14:>15.3f}{15:>15.3f}{16:>15.3f}{17:>11.0f}{9:>15.3f}{10:>15.3f}{11:>15.3f}{18:>15.6f}{19:>15.6f}{14:>15.3f}{15:>15.3f}{16:>15.3f}{17:>11.0f}'.format(
                j+1, t[j],dt[j]/(24*3600.),sid[j],wa_num,1,0,0,0,0,0,0,cos( feed_angle[j]),sin( feed_angle[j]),0,0,1,refant,cos(-feed_angle[j]),sin(-feed_angle[j])) 
                else: 
                    print >> w, '{0:8.0f}{1:>24s}{1:>15s}{1:>11s}{1:>11s}{1:>11s}{1:>11s}{1:>15s}{1:>11s}{1:>15s}{1:>15s}{1:>15s}{2:>15.6f}{3:>15.6f}{4:>15.3f}{5:>15.3f}{6:>15.3f}{7:>11.0f}{1:>15s}{1:>15s}{1:>15s}{8:>15.6f}{9:>15.6f}{4:>15.3f}{5:>15.3f}{6:>15.3f}{7:>11.0f}'.format(
                j+1,"''",cos( feed_angle[j]),sin( feed_angle[j]),0,0,1,refant,cos(-feed_angle[j]),sin(-feed_angle[j]))
        print >> w, '***END*PASS***'
    
    ## then load in with TBIN
    runtbin(indata,'wa_pang.TBIN') 
    # and apply solution SN_OUT + CL_IN = CL_OUT
    # needs to be new SN table 
    runclcal(indata,sn_out,cl_in,cl_out,refant) 
    return 1

################################################

if __name__=='__main__':
    main()

################################################
