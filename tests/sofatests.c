#include "sofa.h"
#include "stdio.h"
#include "stdlib.h"
#include "strings.h"
#define PI 3.1415926535897931

//Helpers
void printmat(char* title, double mat[3][3]) {
    printf("%s:\n[[ %e,%e,%e]\n [%e,%e,%e]\n[%e,%e,%e]]\n",title,
            mat[0][0],mat[0][1],mat[0][2],
            mat[1][0],mat[1][1],mat[1][2],
            mat[2][0],mat[2][1],mat[2][2]);
}

void printvec(char* title, double vec[3]) {
    printf("%s:\n[%e,%e,%e]\n",title,
            vec[0],vec[1],vec[2]);
}

void deg2dms(double deg,int* d,int* m,double* s) {
    *d = (int)(deg);
    *m = (int)((deg-*d)*60);
    *s = ((deg-*d-*m/60.)*3600);
}

//Test functions
void trans_matricies(double ep, double jd1) {
    double dpsi,deps,epsa,rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3],rc2i[3][3];

    iauNut00b(jd1,0,&dpsi,&deps);
    iauPn00(jd1,0,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    printmat("B",rb);
    printmat("P",rp);
    printmat("N",rn);
    printmat("BPN",rbpn);

    iauC2i00b(jd1,0,rc2i);
    printmat("C",rc2i);

}

void trans_coords(double ep, double jd) {
    double vec[3],transvec[3],rc2i[3][3],rbpn[3][3],rc2t[3][3];
    double theta,phi,rasec,decsec;
    int rahr,ramin,decdegi,decmin;
    
    double lng = 10; //ra
    double lat = 20; //dec
    
    iauS2c(lng*PI/180,lat*PI/180,vec);
    iauC2i00b(jd,0,rc2i);
    //printmat("C2I",rc2i);
    //printvec("cartesian vec",vec);
    iauPnm00b(jd, 0, rbpn);
    //printmat("BPN",rbpn);
    //printvec("cartesian vec",vec);
    iauC2t00b(jd,0,jd,0,0,0,rc2t);
    //printmat("C2T",rc2t);
    //printvec("cartesian vec",vec);
    
    iauRxp(rc2i,vec,transvec);
    iauC2s(transvec,&theta,&phi);
    
    
    deg2dms(theta*180/PI/15,&rahr,&ramin,&rasec);
    deg2dms(phi*180/PI,&decdegi,&decmin,&decsec);
    printf("CIRS coords: ra:%ih%im%gs dec:%id%im%gs\n",rahr,ramin,rasec,decdegi,decmin,decsec);
    
    iauRxp(rbpn,vec,transvec);
    iauC2s(transvec,&theta,&phi);
    
    deg2dms(theta*180/PI/15,&rahr,&ramin,&rasec);
    deg2dms(phi*180/PI,&decdegi,&decmin,&decsec);
    printf("Equinox/BPN coords: ra:%ih%im%gs dec:%id%im%gs\n",rahr,ramin,rasec,decdegi,decmin,decsec);
    
    iauRxp(rc2t,vec,transvec);
    iauC2s(transvec,&theta,&phi);
    printf("ITRS coords: lat:%.12g long:%.12g\n",phi*180/PI,theta*180/PI);

}

void earth_rotation(double ep, double jd) {
    printf("ERA:%.12g\n",iauEra00(jd,0));
    printf("GAST:%.12g\n",iauGst00b(jd,0));
    printf("GMST:%.12g\n",iauGmst00(jd,0,jd,0));
}

void earth_pv(double ep, double jd) {
    int res;
    double pvh[2][3],pvb[2][3];
    
    res = iauEpv00(jd,0,pvh,pvb);
    
    if (res!=0) {
        printf("JD not in range 1900-2100 CE as expected by SOFA function!");
    }
    
    printf("Heliocentric pos: %g,%g,%g\n",pvh[0][0],pvh[0][1],pvh[0][2]);
    printf("Heliocentric vel: %g,%g,%g\n",pvh[1][0],pvh[1][1],pvh[1][2]);
    
    printf("SS Barycentric pos: %g,%g,%g\n",pvb[0][0],pvb[0][1],pvb[0][2]);
    printf("SS Barycentric vel: %g,%g,%g\n",pvb[1][0],pvb[1][1],pvb[1][2]);
}

int main(int argc, const char* argv[] ){
    if (argc != 4){
        printf("Need command line arguments: testname epoch jd\n");
    } else {
        double epoch,jd;
        char* testname = argv[1];
        int doall = strcmp(testname,"all")==0;
        
        epoch = atof(argv[2]);
        jd = atof(argv[3]);
        
        if (doall|(strcmp(testname,"trans_matricies")==0)) {
            trans_matricies(epoch,jd);
        }   
        if (doall|(strcmp(testname,"trans_coords")==0)) {
            trans_coords(epoch,jd);
        }
        if (doall|(strcmp(testname,"earth_rotation")==0)) {
            earth_rotation(epoch,jd);
        }    
        if (doall|(strcmp(testname,"earth_pv")==0)) {
            earth_pv(epoch,jd);
        } 
    }
    
    return 0;
    
}
