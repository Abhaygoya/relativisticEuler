#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double length = 1.;
double time = .4;
double gam = 1.67;
int N;
double dx, dt, currentTime;


struct vector {
	double mass;
	double mom;
	double momY;
	double ene;
};

struct primitive {
	double rho;
	double v;
	double eps;
	double p;
	double W;
	double vy;
};

double getx(int i){
	double x = length*(i%N+.5)/N;
	return x-.5;
}

double gety(int i){
	double y = length*(i/N+.5)/N;
	return y-.5;
}

double getcs(int i, struct primitive * P){
	double cs = sqrt((gam-1.)*gam*P[i].eps/(1.+gam*P[i].eps));
	return cs;
}

double getW(int i, struct primitive * P){
	double W = 1./sqrt(1.-pow(P[i].v,2)-pow(P[i].vy,2));
	return W;
}

double geth(int i, struct primitive * P){
	double h = 1.+P[i].eps+P[i].p/P[i].rho;
	return h;
}

void initialize(struct vector * U, struct primitive * P, int i){
	//Get x value, determine U based on x. All calculated from initial vars rho, p, v, and gam
	double pos = getx(i);
	if(pos<0.){
		P[i].rho = 1.0;
		P[i].v = 0.0;
		P[i].vy = 0.0;
		P[i].W=getW(i,P);
		P[i].p = 1000.0;
		P[i].eps = P[i].p/((gam-1.)*P[i].rho);
		
	}
	else{
		P[i].rho = 1.0;
		P[i].v = 0.0;
		P[i].vy = 0.99;
		P[i].W=getW(i,P);
		P[i].p = 0.01;
		P[i].eps = P[i].p/((gam-1.)*P[i].rho);
	}
	U[i].mass=P[i].rho*P[i].W;
	double h = geth(i,P);
	U[i].mom=P[i].rho*h*pow(P[i].W,2)*P[i].v;
	U[i].momY=P[i].rho*h*pow(P[i].W,2)*P[i].vy;
	U[i].ene=P[i].rho*h*pow(P[i].W,2)-P[i].p-U[i].mass;

}

void updateFx(struct vector * U, struct vector * F, struct primitive * P){
	int i;
	for(i=0;i<N*N;i++){
		F[i].mass = U[i].mass*P[i].v;
		F[i].mom = U[i].mom*P[i].v+P[i].p;
		F[i].momY = U[i].momY*P[i].v;
		F[i].ene = U[i].mom-U[i].mass*P[i].v;
	}
}

void updateFy(struct vector * U, struct vector * F, struct primitive * P){
	int i;
	for(i=0;i<N*N;i++){
		F[i].mass = U[i].mass*P[i].vy;
		F[i].mom = U[i].mom*P[i].vy;
		F[i].momY = U[i].momY*P[i].vy+P[i].p;
		F[i].ene = U[i].mom-U[i].mass*P[i].v;
	}
}

double max(double a, double b){
	if(a>=b) return a;
	else return b;
}

double min(double a, double b){
	if(a<=b) return a;
	else return b;
}
//Calculate alphas at i, plus parameter determines alpha+ for 1 and alpha- for 0
double get_alphaX(int i, struct vector * U, struct vector * F, struct primitive * P, int plus){
	double v = .5*(P[i].v+P[i+1].v);
	double cs = .5*(getcs(i,P)+getcs(i+1,P));

	double alpha;
	if(plus==1) alpha = max(0,((v+cs)/(1.+v*cs)));
	else alpha = min(0,((v-cs)/(1.-v*cs)));
	
	//printf("%f", alpha);
	return alpha;	
}

struct vector getFluxX(int i, struct vector * U, struct vector * F, struct primitive * P){
	double alpha_plus = get_alphaX(i,U,F,P,1);
	double alpha_minus = get_alphaX(i,U,F,P,0);
	
	struct vector fhll;
	fhll.mass = (alpha_plus*F[i].mass-alpha_minus*F[i+1].mass+alpha_plus*alpha_minus*(U[i+1].mass-U[i].mass))/(alpha_plus-alpha_minus);
	fhll.mom = (alpha_plus*F[i].mom-alpha_minus*F[i+1].mom+alpha_plus*alpha_minus*(U[i+1].mom-U[i].mom))/(alpha_plus-alpha_minus);
	fhll.momY = (alpha_plus*F[i].momY-alpha_minus*F[i+1].momY+alpha_plus*alpha_minus*(U[i+1].momY-U[i].momY))/(alpha_plus-alpha_minus);
	fhll.ene = (alpha_plus*F[i].ene-alpha_minus*F[i+1].ene+alpha_plus*alpha_minus*(U[i+1].ene-U[i].ene))/(alpha_plus-alpha_minus);
	
	return fhll;
	
}

double get_alphaY(int i, struct vector * U, struct vector * F, struct primitive * P, int plus){
	double v = .5*(P[i].v+P[i+N].v);
	double cs = .5*(getcs(i,P)+getcs(i+N,P));

	double alpha;
	if(plus==1) alpha = max(0,((v+cs)/(1.+v*cs)));
	else alpha = min(0,((v-cs)/(1.-v*cs)));
	
	//printf("%f", alpha);
	return alpha;	
}

struct vector getFluxY(int i, struct vector * U, struct vector * F, struct primitive * P){
	double alpha_plus = get_alphaY(i,U,F,P,1);
	double alpha_minus = get_alphaY(i,U,F,P,0);
	
	struct vector fhll;
	fhll.mass = (alpha_plus*F[i].mass-alpha_minus*F[i+N].mass+alpha_plus*alpha_minus*(U[i+N].mass-U[i].mass))/(alpha_plus-alpha_minus);
	fhll.mom = (alpha_plus*F[i].mom-alpha_minus*F[i+N].mom+alpha_plus*alpha_minus*(U[i+N].mom-U[i].mom))/(alpha_plus-alpha_minus);
	fhll.momY = (alpha_plus*F[i].momY-alpha_minus*F[i+N].momY+alpha_plus*alpha_minus*(U[i+N].momY-U[i].momY))/(alpha_plus-alpha_minus);
	fhll.ene = (alpha_plus*F[i].ene-alpha_minus*F[i+N].ene+alpha_plus*alpha_minus*(U[i+N].ene-U[i].ene))/(alpha_plus-alpha_minus);
	
	return fhll;
	
}

void updatePrimitive(struct vector * U, struct primitive * P){
	int i;
	for(i=0;i<N*N;i++){
		double delta = 10;
		P[i].p=10.0;
		while(fabs(delta)>.00000001){
			//Sets a minimum value for P so NR doesn't choose an impossible guess
			if(P[i].p<(fabs(U[i].mom)-U[i].ene-U[i].mass)) P[i].p=fabs(U[i].mom)-U[i].ene-U[i].mass;
			P[i].v=U[i].mom/(U[i].ene+U[i].mass+P[i].p);
			P[i].vy=U[i].momY/(U[i].ene+U[i].mass+P[i].p);
			P[i].W=getW(i,P);

			P[i].rho=U[i].mass/P[i].W;
			P[i].eps=(U[i].ene+U[i].mass*(1.-P[i].W)+P[i].p*(1.-pow(P[i].W,2)))/(U[i].mass*P[i].W);
			//P(EOS)-P(GUESS) = 0 is the basic f we use; EOS: P = (gam-1)*rho*eps
			double f = (gam-1.)*P[i].rho*P[i].eps-P[i].p;
			double cs = getcs(i,P);
			double dfdx = pow(P[i].v,2)*pow(cs,2)-1.;
			delta = -f/dfdx;
			P[i].p+=delta;
		}
		//To re-update other prim vars after final p is found
		P[i].v=U[i].mom/(U[i].ene+U[i].mass+P[i].p);
		P[i].vy=U[i].momY/(U[i].ene+U[i].mass+P[i].p);
		P[i].W=getW(i,P);
		P[i].rho=U[i].mass/P[i].W;
		P[i].eps=(U[i].ene+U[i].mass*(1.-P[i].W)+P[i].p*(1.+pow(P[i].W,2)))/(U[i].mass*P[i].W);
	
	
	}
}
//Lock the right/left as IC and use 0 gradient on top/bottom
void fixBC(struct vector * U, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		initialize(U,P,i*N+0);
		initialize(U,P,i*N+N-1);
	}
	
	for(i=0;i<N;i++){
		U[i]=U[i+N];
		P[i]=P[i+N];
		U[(N-1)*N+i]=U[(N-2)*N+i];
		P[(N-1)*N+i]=P[(N-2)*N+i];
	}
}
//For second test, Reimann problem in Y direction
void fixBCy(struct vector * U, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		initialize(U,P,i);
		initialize(U,P,i+N*(N-1));
	}
	
	for(i=0;i<N;i++){
		U[i*N]=U[i*N+1];
		P[i*N]=P[i*N+1];
		U[i*N+N-1]=U[i*N+N-2];
		P[i*N+N-1]=P[i*N+N-2];
		P[(N-1)*N+i]=P[(N-2)*N+i];
	}
}

//Fiph arrays contain fluxes through interfaces. Only need N*(N-1) elements technically, but N*N makes indexing easy.
void advanceTime(struct vector * U, struct vector * Fx, struct vector * Fy, struct primitive * P){
	struct vector FiphX[N*N];
	struct vector FiphY[N*N];
	int i;
	int j;
	
	for(j=0;j<N;j++) for(i=0;i<N-1;i++) FiphX[j*N+i] = getFluxX(j*N+i,U,Fx,P);
	for(j=0;j<N-1;j++) for(i=0;i<N;i++) FiphY[j*N+i] = getFluxY(j*N+i,U,Fy,P);
	for(j=0;j<N;j++){
		for(i=0;i<N-1;i++){
			U[j*N+i].mass -= dt*FiphX[j*N+i].mass/dx;
			U[j*N+i].mom -= dt*FiphX[j*N+i].mom/dx;
			U[j*N+i].momY -= dt*FiphX[j*N+i].momY/dx;
			U[j*N+i].ene -= dt*FiphX[j*N+i].ene/dx;
			U[j*N+i+1].mass += dt*FiphX[j*N+i].mass/dx;
			U[j*N+i+1].mom += dt*FiphX[j*N+i].mom/dx;
			U[j*N+i+1].momY += dt*FiphX[j*N+i].momY/dx;
			U[j*N+i+1].ene += dt*FiphX[j*N+i].ene/dx;
		}
	}
	
	for(j=0;j<N-1;j++){
		for(i=0;i<N;i++){
			U[j*N+i].mass -= dt*FiphY[j*N+i].mass/dx;
			U[j*N+i].mom -= dt*FiphY[j*N+i].mom/dx;
			U[j*N+i].momY -= dt*FiphY[j*N+i].momY/dx;
			U[j*N+i].ene -= dt*FiphY[j*N+i].ene/dx;
			U[j*N+i+N].mass += dt*FiphY[j*N+i].mass/dx;
			U[j*N+i+N].mom += dt*FiphY[j*N+i].mom/dx;
			U[j*N+i+N].momY += dt*FiphY[j*N+i].momY/dx;
			U[j*N+i+N].ene += dt*FiphY[j*N+i].ene/dx;
		}
	}
	
	updatePrimitive(U,P); 
	fixBC(U,P);
	//fixBCy(U,P);
	
	updateFx(U,Fx,P);
	updateFy(U,Fy,P);

}

		


	
int main(){
	//Loop through different N's, advancing to time t=.25, and output
	N = 200;
	while(N<201){
		dt = time/(8*N);
		dx = length/N;
	
		struct vector U[N*N];
		struct vector Fx[N*N];
		struct vector Fy[N*N];
		struct primitive P[N*N];

		int i;
		for(i=0;i<N*N;i++) initialize(U,P,i);
		
		updateFx(U,Fx,P);
		updateFy(U,Fy,P);
		
		currentTime=0.;
		while(currentTime<time){
			advanceTime(U,Fx,Fy,P);
			currentTime += dt;
			if(isnan(P[250].rho)==1){
				printf("%f\n",currentTime);
				break;
			}
		}
		
		FILE* f;
		char string[15];
		sprintf(string,"%d",N);
		strcat(string,"_2d.txt");
		f = fopen(string,"w");
		
		for(i=0;i<N*N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f    %15.15f\n", getx(i),gety(i),P[i].rho,P[i].v,P[i].p,P[i].vy);
		//for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),P[i].rho,P[i].W,P[i].p,P[i].v);
		
		
		
		
		fclose(f);
		N = 2*N;
		
	}
	return 0;
}
