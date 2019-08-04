/*

	cube-viewer
	
	To view a 3D cube constructed out of a stack
	of quadratic 24 bit bitmaps
	
	Observer position can be varied
	Shading by dimming by distance to observer
	and dimming by normal vector on 2x2 grid
	pixels
	
	Marc Meidlinger, August 2019
	
*/

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdint.h"


// constants

enum {
	CMD_NOCOMMAND=0,CMD_CONSTRUCT,CMD_VIEW
};

enum {
	OBSALLE=-1,CORNER1=1,CORNER2,CORNER3,
	CORNER4,CORNER5,CORNER6,CORNER7,CORNER8,
	AXXPLUS,AXXMINUS,AXYPLUS,AXYMINUS,AXZPLUS,
	AXZMINUS,CORNERPERS
};

const int32_t OBSMIN=CORNER1;
const int32_t OBSMAX=AXZMINUS;


// classes

struct Bitmap {
	int32_t xlen,ylen,bytes,ybytes;
	unsigned char* bmp;

	void disp(void);

	Bitmap();
	virtual ~Bitmap();
	int32_t setlenxy(const int32_t,const int32_t);
	void save(const char*);
	void load(const char*);
	void setpixel(const int32_t,const int32_t,const int32_t,const int32_t,const int32_t);
	void getpixel(const int32_t,const int32_t,int32_t& r,int32_t& g,int32_t& b);
	void fill(const int32_t,const int32_t,const int32_t);
};

struct Coord {
	float x,y,z;
	int32_t intersectschritt; 
		
	Coord(const Coord&);
	Coord(const double,const double,const double);
	Coord();
	double norm(void);
	double normQ(void);
	char* strI(char*);
	friend Coord operator+(Coord lhs,const Coord& rhs) { lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; return lhs;	}
	friend Coord operator-(Coord lhs,const Coord& rhs) { lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; return lhs;	}
	Coord& operator=(const Coord&);
	void normiere(void);
	void mult(const double);
};

struct RGB {
	unsigned char B,G,R;
	
	void set(const int32_t,const int32_t,const int32_t);
};

typedef RGB *PRGB;

struct BitcubeFull {
	int32_t boxlenx,boxleny,boxlenz; // die aktuelle im Speicher befindliche Box in ihrer Größe
	int32_t boxx0,boxx1,boxy0,boxy1,boxz0,boxz1;
	int32_t anzrgbsproplane;
	Coord boxCubeM;
	double boxcubed;
	PRGB* rgbBmpByZ;
		
	BitcubeFull();
	virtual ~BitcubeFull();
	void disp(void);
	void setlen(const int32_t,const int32_t,const int32_t);
	void fill(const int32_t,const int32_t,const int32_t);
	void save(const char*);
	int32_t load(const char*);
	void getRgbMem(void);
	void setBmpAtZ(Bitmap&,const int32_t);
	void setpixel(const int32_t,const int32_t,const int32_t,const int32_t,const int32_t,const int32_t);
	void getpixel(const int32_t,const int32_t,const int32_t,int32_t&,int32_t&,int32_t&);
	int32_t within(Coord&,int32_t&,int32_t&,int32_t&);
	int32_t within(Coord&);
	int32_t gibtsschon(const int32_t,const int32_t,const int32_t);
};

struct Screen {
	Bitmap bmp;
	int32_t lenx,leny;
	double abstandDelta;
	double normVektor;
	Coord* coords; 
	Coord vektor;
	Coord scrxv,scryv;
	RGB transp0,transp1; 
	double obsabst;
	Coord observer;
		
	Screen();
	virtual ~Screen();
	void setScreenSize(const int32_t,const int32_t);
	void setObserver(Coord&,Coord&,Coord&,Coord&);
	void slide(BitcubeFull*);
	int32_t getIdx(const int32_t,const int32_t);
	void setTransparent(RGB&,RGB&);
	void dimmNyDistance(const double dimmfaktor,const int32_t s0,const int32_t s1);
	void dimmByNormal(void);
	double distanceToScreen(Coord&);
};


// forward declarations

void swapCoord(Coord&,Coord&);
void crossproduct(Coord&,Coord&,Coord&);
double scalarproduct(Coord&,Coord&);
char* upper(char*);
inline int32_t maximumI(const int32_t,const int32_t);
inline int32_t minimumI(const int32_t,const int32_t);
inline double minimumD(const double,const double);


// globals

char ccbfn[1024];
int32_t obspos0,obspos1;
Coord userobserver;
BitcubeFull bitcube;
RGB transpcolor0,transpcolor1;
Screen screen;
Coord obsUL,obx,oby;
double mindistancetoscreen;


// struct bitmap

void Bitmap::load(const char* fn) {
	disp();
	FILE *fbmp=fopen(fn,"rb");
	if (!fbmp) return;
	
	int32_t w;
	fseek(fbmp,18,SEEK_SET);
	int32_t xl,yl;
	fread(&w,sizeof(w),1,fbmp); xl=w;
	fread(&w,sizeof(w),1,fbmp); yl=w;
	setlenxy(xl,yl);
	fseek(fbmp,54,SEEK_SET);
	fread(bmp,sizeof(unsigned char),bytes,fbmp);

	fclose(fbmp);
}

void Bitmap::setpixel(const int32_t x,const int32_t y,const int32_t r,const int32_t g,const int32_t b) {
	int32_t offset=y*ybytes+3*x;
	if ((offset<0) || ((offset+2)>=bytes)) return;
	bmp[offset]=b;
	bmp[offset+1]=g;
	bmp[offset+2]=r;
}

void Bitmap::getpixel(const int32_t x,const int32_t y,int32_t& r,int32_t& g,int32_t& b) {
	int32_t offset=y*ybytes+3*x;
	if ((offset<0) || ((offset+2)>=bytes)) return;
	b=bmp[offset];
	g=bmp[offset+1];
	r=bmp[offset+2];
}

int32_t Bitmap::setlenxy(const int32_t xl,const int32_t yl) {
	if ((!bmp)||(xlen!=xl)||(ylen!=yl)) {
		disp();
		xlen=xl; ylen=yl;
		ybytes=3*xlen;
		bytes=ylen*ybytes;
		bmp=new unsigned char[bytes];
		if (!bmp) return 0;
	}
	
	return 1;
}

void Bitmap::save(const char* fn) {
	FILE *fbmp=fopen(fn,"wb");
	unsigned char header1[18]={
		0x42,0x4D,0xF6,0xC6,0x2D,0x00,0x00,0x00,
		0x00,0x00,0x36,0x00,0x00,0x00,0x28,0x00,
		0x00,0x00
	};
	fwrite(header1,sizeof(unsigned char),18,fbmp);
	uint32_t w=xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w=ylen; 
	fwrite(&w,sizeof(w),1,fbmp);
	unsigned char header2[]={
		0x01,0x00,0x18,0x00,0x00,0x00,0x00,0x00,
		0xC0,0xC6,0x2D,0x00,0xC4,0x0E,0x00,0x00,
		0xC4,0x0E,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00
	};
	fwrite(header2,sizeof(unsigned char),28,fbmp);
	fwrite(bmp,bytes,1,fbmp);

	fclose(fbmp);
}

Bitmap::Bitmap(void) {
	xlen=ylen=bytes=ybytes=0;
	bmp=NULL;
}

void Bitmap::disp(void) {
	if (bmp) { 
		delete[] bmp; 
		bmp=NULL; 
		bytes=xlen=ylen=0; 
	}
}

Bitmap::~Bitmap(void) {
	disp();
}

void Bitmap::fill(const int32_t r,const int32_t g,const int32_t b) {
	int32_t o;
	for(int32_t y=0;y<ylen;y++) {
		o=y*ybytes;
		for(int32_t x=0;x<xlen;x++) {
			bmp[o]=b; bmp[o+1]=g; bmp[o+2]=r;
			o+=3;
		}
	}
}


// struct rgb

void RGB::set(const int32_t ar,const int32_t ag,const int32_t ab) {
	R=ar; G=ag; B=ab;
}


// struct bitcubefull

void BitcubeFull::fill(const int32_t ar,const int32_t ag,const int32_t ab) {
	for(int32_t z=0;z<boxlenz;z++) {
		if (rgbBmpByZ[z]) {
			for(int32_t i=0;i<anzrgbsproplane;i++) {
				rgbBmpByZ[z][i].R=ar;
				rgbBmpByZ[z][i].G=ag;
				rgbBmpByZ[z][i].B=ab;
			}
		}
	}
}

BitcubeFull::BitcubeFull() {
	anzrgbsproplane=0;
	boxlenx=boxleny=boxlenz=0;
	rgbBmpByZ=NULL;
}

BitcubeFull::~BitcubeFull() {
	disp();
}

void BitcubeFull::disp(void) {
	if (rgbBmpByZ) {
		for(int32_t z=0;z<boxlenz;z++) {
			if (rgbBmpByZ[z]) delete[] rgbBmpByZ[z];
		}
		delete[] rgbBmpByZ;
		rgbBmpByZ=NULL;
		anzrgbsproplane=0;
	}
}

void BitcubeFull::setBmpAtZ(Bitmap& q,const int32_t az) {
	memcpy(&rgbBmpByZ[az][0],q.bmp,q.xlen*q.ylen*3);
}

void BitcubeFull::setlen(const int32_t ax,const int32_t ay,const int32_t az) {
	if (rgbBmpByZ) disp();
	
	boxx1=boxy1=boxz1=-1;
	boxx0=boxy0=boxz0=99999999;
	boxlenx=ax;
	boxleny=ay;
	boxlenz=az;
	anzrgbsproplane=ax*ay;
	boxcubed=sqrt(boxlenx*boxlenx+boxleny*boxleny+boxlenz*boxlenz);
	boxCubeM.x=(boxlenx >> 1);
	boxCubeM.y=(boxleny >> 1);
	boxCubeM.z=(boxlenz >> 1);
	
	getRgbMem();
}

void BitcubeFull::getRgbMem(void) {
	rgbBmpByZ=new PRGB[boxlenz];
	if (!rgbBmpByZ) {
		printf("Memory error.\n");
		exit(99);
	}
	for(int32_t z=0;z<boxlenz;z++) {
		rgbBmpByZ[z]=new RGB[anzrgbsproplane];
		if (!rgbBmpByZ[z]) {
			printf("Memory error.\n");
			exit(99);
		}
	}
	
	anzrgbsproplane=boxlenx*boxleny;
}

int32_t BitcubeFull::load(const char* afn) {
	FILE *f=fopen(afn,"rb");
	if (!f) return -1;
	
	disp();
	fread(&boxlenx,sizeof(int),1,f);
	fread(&boxleny,sizeof(int),1,f);
	fread(&boxlenz,sizeof(int),1,f);
	anzrgbsproplane=boxlenx*boxleny;
	boxcubed=sqrt(boxlenx*boxlenx+boxleny*boxleny+boxlenz*boxlenz);
	boxCubeM.x=(boxlenx >> 1);
	boxCubeM.y=(boxleny >> 1);
	boxCubeM.z=(boxlenz >> 1);
	getRgbMem();
	for(int32_t z=0;z<boxlenz;z++) {
		fread(&rgbBmpByZ[z][0],sizeof(RGB),anzrgbsproplane,f);
	}
	
	fclose(f);

	// where do the actual points lie ?
	boxx0=boxlenx;
	boxy0=boxleny;
	boxz0=boxlenz;
	boxx1=boxy1=boxz1=-1;
	for(int32_t z=0;z<boxlenz;z++) for(int32_t y=0;y<boxleny;y++) for(int32_t x=0;x<boxlenx;x++) {
		if (rgbBmpByZ[z][y*boxlenx+x].R > 1) {
			if (x<boxx0) boxx0=x;
			if (y<boxy0) boxy0=y;
			if (z<boxz0) boxz0=z;
			if (x>boxx1) boxx1=x;
			if (y>boxy1) boxy1=y;
			if (z>boxz1) boxz1=z;
		}
	}

	if (boxx0 >= boxlenx) boxx0=0;
	if (boxy0 >= boxleny) boxy0=0;
	if (boxz0 >= boxlenz) boxz0=0;
	if (boxx1 < 0) boxx1=boxlenx;
	if (boxy1 < 0) boxy1=boxleny;
	if (boxz1 < 0) boxz1=boxlenz;
	
	return 1;
}

void BitcubeFull::save(const char* fn) {
	FILE *f=fopen(fn,"wb");
	if (!f) return;
	
	fwrite(&boxlenx,sizeof(int),1,f);
	fwrite(&boxleny,sizeof(int),1,f);
	fwrite(&boxlenz,sizeof(int),1,f);
	for(int32_t z=0;z<boxlenz;z++) {
		fwrite(&rgbBmpByZ[z][0],sizeof(RGB),anzrgbsproplane,f);
	}
	
	fclose(f);
}

void BitcubeFull::getpixel(const int32_t px,const int32_t py,const int32_t pz,int32_t& r,int32_t& g,int32_t& b) {
	int32_t off=py*boxlenx+px;
	r=rgbBmpByZ[pz][off].R;
	g=rgbBmpByZ[pz][off].G;
	b=rgbBmpByZ[pz][off].B;
}

void BitcubeFull::setpixel(const int32_t px,const int32_t py,const int32_t pz,const int32_t ar,const int32_t ag,const int32_t ab) {
	int32_t off=py*boxlenx+px;
	rgbBmpByZ[pz][off].set(ar,ag,ab);
	if (px < boxx0) boxx0=px;
	if (py < boxy0) boxy0=py;
	if (pz < boxz0) boxz0=pz;
	if (px > boxx1) boxx1=px;
	if (py > boxy1) boxy1=py;
	if (pz > boxz1) boxz1=pz;
}

int32_t BitcubeFull::within(Coord& ap,int32_t& ar,int32_t& ag,int32_t& ab) {
	if (
		(ap.x>=0) && (ap.x<boxlenx)
		&& (ap.y>=0) && (ap.y<boxleny)
		&& (ap.z>=0) && (ap.z<boxlenz)
	) {
		int32_t idx=(int)(floor(ap.y)*boxlenx+floor(ap.x));
		int32_t tmpz=(int)floor(ap.z);
		
		ar=rgbBmpByZ[tmpz][idx].R;
		ag=rgbBmpByZ[tmpz][idx].G;
		ab=rgbBmpByZ[tmpz][idx].B;

		return 1;
	}
	
	return 0;
}

int32_t BitcubeFull::within(Coord& ap) {
	// no color necessary
	if (
		(ap.x>=0) && (ap.x<boxlenx)
		&& (ap.y>=0) && (ap.y<boxleny)
		&& (ap.z>=0) && (ap.z<boxlenz)
	) return 1;
	
	return 0;
}

// struct coord

double Coord::norm(void) { 
	return sqrt(x*x+y*y+z*z); 
}

double Coord::normQ(void) { 
	return (x*x+y*y+z*z); 
}

Coord& Coord::operator=(const Coord& a) {
	if (this != &a) {
		x=a.x;
		y=a.y;
		z=a.z;
	}
	return *this;
}
 
void Coord::normiere(void) {
	double d=1.0; 
	d /= norm();
	mult(d);
}

void Coord::mult(const double d) {
	x *= d; y *= d;	z *= d;
}

Coord::Coord() {
	x=y=z=0.0;
}

Coord::Coord(const double ax,const double ay,const double az) {
	x=ax; y=ay;	z=az;
}

Coord::Coord(const Coord& a) {
	x=a.x; y=a.y; z=a.z;
}

char* Coord::strI(char* erg) {
	sprintf(erg,"(%.0lf,%.0lf,%.0lf)",x,y,z);
	return erg;
}

// struct screen

double Screen::distanceToScreen(Coord& p) {
	return ( (vektor.x*p.x + vektor.y*p.y + vektor.z*p.z - abstandDelta) / normVektor);
}

Screen::Screen() {
	coords=NULL;
	lenx=leny=0;
}

Screen::~Screen() {
	if (coords) { 
		delete[] coords;
		coords=NULL;
		lenx=leny=0;
	}
}

void Screen::setTransparent(RGB& a0,RGB& a1) {
	transp0.R=a0.R;
	transp0.G=a0.G;
	transp0.B=a0.B;
	transp1.R=a1.R;
	transp1.G=a1.G;
	transp1.B=a1.B;
}

int32_t Screen::getIdx(const int32_t x,const int32_t y) {
	return (y*bmp.xlen+x);
}

void Screen::setScreenSize(const int32_t axl,const int32_t ayl) {
	if (coords) delete[] coords;
	
	coords=new Coord[axl*ayl];
	if (!coords) {
		printf("Memory error.\n");
		exit(99);
	}

	lenx=axl;
	leny=ayl;
	bmp.setlenxy(axl,ayl);
}

void Screen::setObserver(Coord& obsUL,Coord& rvx,Coord& rvy,Coord& Mcube) {
	if (!coords) return;

	rvx.normiere();
	rvy.normiere();
	scrxv=rvx;
	scryv=rvy;
	
	const double x0=0.5*bmp.xlen;
	const double y0=0.5*bmp.ylen;
	
	vektor.x=(Mcube.x-(obsUL.x+x0*rvx.x+y0*rvy.x));
	vektor.y=(Mcube.y-(obsUL.y+x0*rvx.y+y0*rvy.y));
	vektor.z=(Mcube.z-(obsUL.z+x0*rvx.z+y0*rvy.z));
	vektor.normiere();
	
	// where does the screen lie currently in
	// the cube ?
	int32_t offset=0;
	for(int32_t y=0;y<bmp.ylen;y++) for(int32_t x=0;x<bmp.xlen;x++) {
		coords[offset].x=obsUL.x+x*rvx.x+y*rvy.x;
		coords[offset].y=obsUL.y+x*rvx.y+y*rvy.y;
		coords[offset].z=obsUL.z+x*rvx.z+y*rvy.z;
		coords[offset].intersectschritt=-1; 
		offset++;
	}
	
	abstandDelta=vektor.x*obsUL.x + vektor.y*obsUL.y + vektor.z*obsUL.z;
	normVektor=vektor.norm();
}

void Screen::dimmNyDistance(const double dimmfaktor,const int32_t s0,const int32_t s1) {
	const int32_t XLEN=bmp.xlen;
	const int32_t YLEN=bmp.ylen;

	int32_t breite=s0-s1; 
	double df=1.0; 
	if (breite<=0) breite=1;
	
	df /= breite;
	df *= dimmfaktor;

	int32_t offset=0;

	for(int32_t y=0;y<YLEN;y++) {
		for(int32_t x=0;x<XLEN;x++) {
			if (coords[offset].intersectschritt>0) { 
				int32_t r=0,g=0,b=0;
				bmp.getpixel(x,y,r,g,b);
				double dd=1-df*(s0 - coords[offset].intersectschritt);
				bmp.setpixel(x,y,(int)floor(dd*r),(int)floor(dd*g),(int)floor(dd*b));
			} // dimm point
			
			offset++;
		} // x
	} // y
}

// principal function

void Screen::slide(BitcubeFull* pbc) {
	if (!pbc) return;
	
	double len=sqrt(pbc->boxlenx*pbc->boxlenx+pbc->boxleny*pbc->boxleny+pbc->boxlenz*pbc->boxlenz);
	int32_t schritte=2*(int)len;
	int32_t MAXSCR=(int)ceil(1.2*sqrt(pbc->boxlenx*pbc->boxlenx+pbc->boxleny*pbc->boxleny+pbc->boxlenz*pbc->boxlenz)); // diese maximale Richtung kann irgendein POunkt, der den
	
#define CUBAUS(PTR) \
( (PTR->x<0) || (PTR->y<0) || (PTR->z<0) || (PTR->x>=pbc->lenx) || (PTR->y>=pbc->leny) || (PTR->z>=pbc->lenz) )		

	int32_t maxanzschritte=2*(int)len;
	schritte=maxanzschritte;
	int32_t intersectidx0=-1,intersectidx1=-1; // invers, d.h. je kleiner Index, desto näher an Cube
	int32_t offset=-1;
	int32_t erster=1;
	const int32_t XLEN=bmp.xlen;
	const int32_t YLEN=bmp.ylen;
	
	const int32_t NOCH0=64;
	int32_t noch=NOCH0;
	
	// clear screen from previous intersections if any
	offset=0;
	for(int32_t y=0;y<YLEN;y++) {
		for(int32_t x=0;x<XLEN;x++) {
			coords[offset].intersectschritt=-1;
			offset++;
		}
	}
	
	int32_t sprung=(int)floor(0.9*mindistancetoscreen);

	const int32_t SPRUNGSCHRITT=8;
	const int32_t ZURUECKSCHRITTE=64;
	
	int32_t inkubus=1;
	Coord add;
	
	// one big jump, maybe alread into the points ?
	while (inkubus>0) {
		inkubus=0;
		offset=0;
		add=vektor;
		add.mult(sprung);
		for(int32_t y=0;y<YLEN;y++) {
			for(int32_t x=0;x<XLEN;x++) {
				if (coords[offset].intersectschritt < 0) { 
					Coord temp=coords[offset];
					temp=temp+add;
					if (pbc->within(temp) > 0) {
						// step a bit back
						sprung-=ZURUECKSCHRITTE;
						inkubus=1;
						break;
					}
				}
				offset++;
			} // for x
			
			if (inkubus>0) break;
		} // for y
	} // while
	
	// now in smaller jumps to get closer without
	// having to check every step
	int32_t zaehle=schritte;
	while (zaehle>0) {
		int32_t treffer=0;
		
		int32_t offset=0;
		add=vektor; 
		add.mult(sprung);
		for(int32_t y=0;(y<YLEN);y++) {
			for(int32_t x=0;(x<XLEN);x++) {
				if (coords[offset].intersectschritt < 0) { 
					Coord temp=coords[offset];
					temp=temp+ add;
					if (pbc->within(temp) > 0) {
						treffer=1;
						break;
					}
				}
				offset++;
			}
			if (treffer>0) break;
		} // fors

		if (treffer>0) {
			// step a bit back
			sprung-=SPRUNGSCHRITT;
			break;
		}
		
		sprung+=SPRUNGSCHRITT;
		zaehle-=SPRUNGSCHRITT;
	}
	
	if (sprung > 0) {
		// adjust alls creen coordinates to the
		// starting position where now every single
		// step walking through space is taken
		offset=0;
		Coord vsprung=vektor;
		vsprung.mult(sprung);
		for(int32_t y=0;y<YLEN;y++) {
			for(int32_t x=0;x<XLEN;x++) {
				coords[offset]=coords[offset]+vsprung;
				offset++;
			}
		}
	}

	schritte-=sprung; 
	
	int32_t rc,gc,bc;
	// screen region to test and adjust coordinates
	int32_t x0=0,x1=XLEN-1,y0=0,y1=YLEN-1;
	
	while (schritte>0) {
		if ( (noch--)== 0) {
			printf("%i ",schritte);
			noch=NOCH0;
		}
		schritte--;
		int32_t scrschneidet=0;
		
		for(int32_t y=y0;y<=y1;y++) {
			offset=y*XLEN + x0;
			for(int32_t x=x0;x<=x1;x++) {
				if (coords[offset].intersectschritt < 0) { 
					if (pbc->within(coords[offset],rc,gc,bc) > 0) {
						// screen intersects with the cube itself
						// not necessarily with a colored pixel though
						scrschneidet=1;
						x0=maximumI(x0,x-MAXSCR);
						x1=minimumI(x1,x+MAXSCR);
						y0=maximumI(y0,y-MAXSCR);
						y1=minimumI(y1,y+MAXSCR);
						
						if (
							( (rc+gc+bc) == 0) ||
							(
								(rc >= transp0.R) && (rc <= transp1.R) &&
								(gc >= transp0.G) && (gc <= transp1.G) &&
								(bc >= transp0.B) && (bc <= transp1.B)
							)
						) {
							coords[offset]=coords[offset]+vektor;
						} else {
							if (intersectidx0<0) intersectidx0=schritte;
							intersectidx1=schritte;
							
							coords[offset].intersectschritt=schritte;
							
							bmp.setpixel(x,y,maximumI(1,rc),gc,bc);
						}
					} else {
						// move that screen pixel to its next location
						coords[offset]=coords[offset]+vektor;
					}
				} 
				
				offset++;
			} // for x

		} // for y
		
		if (scrschneidet<=0) {
			// screen now is behind the whole cube
			if (erster<=0) break; 
		} else if (erster>0) erster=0;
		
	} // while
	
	dimmNyDistance(0.75,intersectidx0,intersectidx1);
	dimmByNormal();
}


// struct screen

void Screen::dimmByNormal(void) {
	// tile the screen into 2x2 pixel grids
	// take all normal vector of those pixel
	// that are colored and then their average
	// angle between that vector and the observer
	// determines dimming
	// only works for 3 or 4 pixels colored,
	// in the other two cases: no dimming here
	
	const double VNORM=vektor.norm();
	
	for(int32_t y=1;y<leny;y++) {
		for(int32_t x=1;x<lenx;x++) {
			int32_t pixele=0;
			int32_t nidx=-1;
			double mindQ=-1;
			Coord ps[4];
			for(int32_t dy=-1;dy<=0;dy++) for(int32_t dx=-1;dx<=0;dx++) {
				const int32_t idx=(y+dy)*lenx+(x+dx);
				if (coords[idx].intersectschritt>=0) {
					ps[pixele]=coords[idx];
					pixele++;
					Coord dist=observer - coords[idx];
					double d=dist.normQ();
					if ((d < mindQ) || (mindQ<-0.5)) {
						mindQ=d;
						nidx=(pixele-1);
					}
				}
			}
			
			if (pixele <= 2) continue;
			
			Coord nvektor;
			if (pixele==3) {
				Coord v1,v2;
				if (nidx != 0) swapCoord(ps[0],ps[nidx]);
				nidx=0;
				v1=ps[1]-ps[0];
				v2=ps[2]-ps[0];
				crossproduct(v1,v2,nvektor);
			} else if (pixele==4) {
				Coord v1,v2,e;
				nvektor=Coord(0,0,0);
				if (nidx != 0) swapCoord(ps[0],ps[nidx]);
				nidx=0;
				
				#define EBENE(AA,BB) \
				v1=ps[AA]-ps[0];\
				v2=ps[BB]-ps[0];\
				crossproduct(v1,v2,e);\
				nvektor=nvektor+e;\
				
				EBENE(1,2)
				EBENE(1,3)
				EBENE(2,3)
			} // 4 pixels
			
			double nvn=nvektor.norm();
			// nil vector ?
			if (fabs(nvn)<1E-5) continue;

			double alpha=
				fabs(
					scalarproduct(vektor,nvektor)
					/ (VNORM*nvn
				)
			);
			alpha=1-0.2*alpha;

			for(int32_t dy=-1;dy<=0;dy++) for(int32_t dx=-1;dx<=0;dx++) {
				const int32_t idx=(y+dy)*lenx+(x+dx);
				if (coords[idx].intersectschritt<0) continue;

				int32_t r=0,g=0,b=0;
				bmp.getpixel(x+dx,y+dy,r,g,b);
				bmp.setpixel(x+dx,y+dy,floor(alpha*r),floor(alpha*g),floor(alpha*b));
			}
		} // x
	} // y
}


// principal functions

void cubing(
	Coord& obsUL,Coord& obx,Coord& oby,
	Screen& scr,
	BitcubeFull& bc,
	RGB& t0,RGB& t1,
	char* fn
) {
	scr.setObserver(obsUL,obx,oby,bc.boxCubeM);
	scr.bmp.fill(0,0,0);
	scr.setTransparent(t0,t1);
	
	mindistancetoscreen=bc.boxcubed*bc.boxcubed;
	
	// distance to all 8 corners of the cube
	for(int32_t z=0;z<=1;z++) 
	for(int32_t y=0;y<=1;y++) 
	for(int32_t x=0;x<=1;x++) {
		Coord temp=Coord(	x*(bc.boxx1-bc.boxx0)+bc.boxx0,
					y*(bc.boxy1-bc.boxy0)+bc.boxy0,
					z*(bc.boxz1-bc.boxz0)+bc.boxz0 );
		mindistancetoscreen=minimumD(
			mindistancetoscreen,
			scr.distanceToScreen(temp));
	}
	
	scr.slide(&bc);
	
	scr.bmp.save(fn);
}

// determine vectors of the screen orientation
void determineScrRV(Screen& scr,Coord& observer,Coord& Pmcube,Coord& obsUL,Coord& rx,Coord& ry) {
	Coord s=Pmcube - observer;
	
	const double fsx=fabs(s.x);
	const double fsy=fabs(s.y);
	const double fsz=fabs(s.z);
	
	if ( (fsx<1E-10) && (fsy<1E-10) && (fsz>1E-10) ) {
		rx=Coord(1,0,0); ry=Coord(0,1,0);
	} else
	if ( (fsx<1E-10) && (fsy>1E-10) && (fsz<1E-10) ) {
		rx=Coord(1,0,0); ry=Coord(0,0,1);
	} else
	if ( (fsx>1E-10) && (fsy<1E-10) && (fsz<1E-10) ) {
		rx=Coord(0,0,1); ry=Coord(0,1,0);
	} else
	if ( (fsx>1E-10) && (fsy>1E-10) && (fsz<1E-10) ) {
		rx=Coord(0,0,1); ry=Coord(1,-s.x/s.y,0);
	} else
	if ( (fsx>1E-10) && (fsy<1E-10) && (fsz>1E-10) ) {
		rx=Coord(0,1,0); ry=Coord(1,0,-s.x/s.z);
	} else
	if ( (fsx<1E-10) && (fsy>1E-10) && (fsz>1E-10) ) {
		rx=Coord(1,0,0); ry=Coord(0,1,-s.y/s.z);
	} else
	if ( (fsx>1E-10) && (fsy>1E-10) && (fsz>1E-10) ) {
		//rx=Coord(1,0,-s.z/s.x); ry=Coord(0,1,-s.y/s.z);
		rx=Coord(-s.y,s.x,0); 
		ry=Coord(
			s.z/(-s.y*s.y/s.x-s.x),
			s.y*s.z/(-s.y*s.y-s.x*s.x),
			1
		);
	} else {
		// nil vector
		printf("Error. Observing nil vector.\n");
		exit(99);
	}
	
	if ( observer.z < -0.01 ) rx.mult(-1); 
	
	rx.normiere(); rx.mult(scr.bmp.xlen);
	ry.normiere(); ry.mult(scr.bmp.ylen);
	
	obsUL.x = observer.x - 0.5*rx.x - 0.5*ry.x;	
	obsUL.y = observer.y - 0.5*rx.y - 0.5*ry.y;	
	obsUL.z = observer.z - 0.5*rx.z - 0.5*ry.z;	
}

void setObserver(BitcubeFull& bc,const int32_t pixel,Coord& observer) {
	switch (pixel) {
		case CORNER1: observer=Coord(-(bc.boxlenx >> 1),bc.boxlenx+(bc.boxleny >> 1),bc.boxlenx+(bc.boxlenz >> 1)); break;
		case CORNER2: observer=Coord(bc.boxlenx+(bc.boxlenx >> 1),bc.boxleny+(bc.boxleny >> 1),bc.boxlenz+(bc.boxlenz >> 1)); break;
		case CORNER3: observer=Coord(bc.boxlenx+(bc.boxlenx >> 1),bc.boxleny+(bc.boxleny >> 1),-(bc.boxlenz >> 1)); break;
		case CORNER4: observer=Coord(-(bc.boxlenx >> 1),bc.boxleny+(bc.boxleny >> 1),-(bc.boxlenz >> 1)); break;
		case CORNER5: observer=Coord(-(bc.boxlenx >> 1),-(bc.boxleny >> 1),bc.boxlenz+(bc.boxlenz >> 1)); break;
		case CORNER6: observer=Coord(bc.boxlenx+(bc.boxlenx >> 1),-(bc.boxleny >> 1),bc.boxlenz+(bc.boxlenz >> 1)); break;
		case CORNER7: observer=Coord(bc.boxlenx+(bc.boxlenx >> 1),-(bc.boxleny >> 1),-(bc.boxlenz >> 1)); break;
		case CORNER8: observer=Coord(-(bc.boxlenx >> 1),-(bc.boxleny >> 1),-(bc.boxlenz >> 1)); break;
		case AXXPLUS: observer=Coord(-(bc.boxlenx >> 1),bc.boxleny >> 1,bc.boxlenz >> 1); break;
		case AXXMINUS: observer=Coord(bc.boxlenx+(bc.boxlenx >> 1),bc.boxleny >> 1,bc.boxlenz >> 1); break;
		case AXYPLUS: observer=Coord(bc.boxlenx >> 1,-(bc.boxleny >> 1),bc.boxlenz >> 1); break;
		case AXYMINUS: observer=Coord(bc.boxlenx >> 1,bc.boxleny+(bc.boxleny >> 1),bc.boxlenz >> 1); break;
		case AXZPLUS: observer=Coord(bc.boxlenx >> 1,bc.boxleny >> 1,-(bc.boxlenz >> 1)); break;
		case AXZMINUS: observer=Coord(bc.boxlenx >> 1,bc.boxleny >> 1,bc.boxlenz+(bc.boxlenz >> 1)); break;
		case CORNERPERS: observer=Coord(0,bc.boxleny,bc.boxlenx+(bc.boxlenz >> 1)); break;
		default: observer=Coord(bc.boxlenx+100,bc.boxlenx+100,bc.boxlenx+100); break;
	}
}


// general routines

inline double minimumD(const double a,const double b) {
	if (a < b) return a;
	return b;
}

inline int32_t maximumI(const int32_t a,const int32_t b) {
	if (a > b) return a;
	
	return b;
}

inline int32_t minimumI(const int32_t a,const int32_t b) {
	if (a < b) return a;
	return b;
}

char* upper(char* s) {
	if (!s) return 0;
	for(uint32_t i=0;i<strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

void swapCoord(Coord& a,Coord& b) {
	Coord c=a;
	a=b;
	b=c;
}

void crossproduct(Coord& a,Coord& b,Coord& e) {
	e.x=a.y*b.z-a.z*b.y;
	e.y=a.z*b.x-a.x*b.z;
	e.z=a.x*b.y-a.y*b.x;
}

double scalarproduct(Coord& a,Coord& b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

int32_t main(int32_t argc, char *argv[]) {
	int32_t cmd=CMD_NOCOMMAND;
	char tmp[1024];
	
	// standard values
	obspos0=OBSMIN;
	obspos1=OBSMAX;
	transpcolor0.set(0,0,0);
	transpcolor1.set(10,10,10);
	strcpy(ccbfn,"bc1.ccb");

	for(int32_t i=1;i<argc;i++) {
		upper(argv[i]);
		int32_t l=strlen(argv[i]);
		if (l>1000) { l=1000; argv[i][1000]=0; }

		if (strstr(argv[i],"CMD=")==argv[i]) {
			if (strstr(&argv[i][4],"CONSTR")==&argv[i][4]) cmd=CMD_CONSTRUCT;
			else cmd=CMD_VIEW;
		} else
		if (strstr(argv[i],"FILE=")==argv[i]) {
			strcpy(ccbfn,&argv[i][5]);
			if (strstr(ccbfn,".CCB") == NULL) sprintf(&ccbfn[strlen(ccbfn)],".CCB");
		} else
		if (strstr(argv[i],"OBSERVER=")==argv[i]) {
			strcpy(tmp,&argv[i][9]);
			// either -1 => all observer position using
			// or 1..14: specific corner/plane points
			// or a coordinate triple
			int32_t x,y,z;
			// standard

			if (sscanf(tmp,"%i,%i,%i",&x,&y,&z) == 3) {
				userobserver=Coord(x,y,z);
				obspos0=obspos1=CORNERPERS;
			} else if (sscanf(tmp,"%i",&x) == 1) {
				if ( (x>=OBSMIN) && (x<=OBSMAX) ) {
					obspos0=obspos1=x;
				} 
			} 
		} else
		if (strstr(argv[i],"TRANSP=")==argv[i]) {
			int32_t r0,g0,b0,r1,g1,b1;
			if (sscanf(
				&argv[i][7],"%i,%i,%i,%i,%i,%i",
				&r0,&g0,&b0,&r1,&g1,&b1
			) == 6) {
				transpcolor0.R=r0; transpcolor0.G=g0; transpcolor0.B=b0;
				transpcolor1.R=r1; transpcolor1.G=g1; transpcolor1.B=b1;
			} 
		}
	} // i
	
	if (cmd==CMD_CONSTRUCT) {
		int32_t lanz=0;
		int32_t geladen=1;
		Bitmap bmp;
		bitcube.fill(0,0,0);
		char tmp[1024];
		
		while (1) {
			printf(".");
			sprintf(tmp,"%08i.bmp",geladen);
			bmp.load(tmp);

			if (lanz==0) {
				if (bmp.xlen != bmp.ylen) {
					printf("Error. Only quadratic images supported.\n");
					exit(99);
				}
				lanz=bmp.xlen;
				bitcube.setlen(bmp.xlen,bmp.ylen,lanz);
			}

			if (bmp.xlen>0)	bitcube.setBmpAtZ(bmp,geladen-1);
		
			geladen++;
			if (geladen > lanz) break;
		}
		bitcube.save(ccbfn);
	} else
	if (cmd==CMD_VIEW) {
		// load bitcube
		if (bitcube.load(ccbfn) <= 0) {
			printf("Error. Butcubefile not found.\n");
			exit(99);
		}
		printf("Bitcube size %i loaded.\n",bitcube.boxlenx);
		
		// slide screen through space
		char fn[1024],posstr[1024];
		Coord observer;
		int32_t filenr=1;
		int32_t c=maximumI(bitcube.boxlenx,maximumI(bitcube.boxleny,bitcube.boxlenz)) << 1;
		screen.setScreenSize(c,c);

		for(int32_t p=obspos0;p<=obspos1;p++) {
			printf("\nobserver #%i ... ",p);
			if (p==CORNERPERS) observer=userobserver; 
			else setObserver(bitcube,p,observer);
				
			// screen moving vector
			determineScrRV(
				screen,
				observer,
				bitcube.boxCubeM,
				obsUL,obx,oby
			);
			
			sprintf(fn,"%s-scr%02i_%s_%04i.bmp",ccbfn,p,observer.strI(posstr),filenr++);
				
			// walk screen through cube to encounter
			// colored points (outside transparency)
			cubing(
				obsUL,obx,oby,
				screen,
				bitcube,
				transpcolor0,transpcolor1,
				fn
			);
			}
		
	} else {
		printf("Unknown command.\n");
		exit(99);
	}
	
	return 0;
}

