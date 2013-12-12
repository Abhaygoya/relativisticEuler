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
	double momT;
	double ene;
};

struct primitive {
	double rho;
	double v;
	double eps;
	double p;
	double W;
	double v_trans;
};

double getx(int i){
	double x = length*(i+.5)/N;
	return x-.5;
}

double getcs(int i, struct primitive * P){
	double cs = sqrt((gam-1.)*gam*P[i].eps/(1.+gam*P[i].eps));
	return cs;
}

double getW(int i, struct primitive * P){
	double W = 1./sqrt(1.-pow(P[i].v,2)-pow(P[i].v_trans,2));
	return W;
}

double geth(int i, struct primitive * P){
	double h = 1.+P[i].eps+P[i].p/P[i].rho;
	return h;
}

void initialize(struct vector * U, struct primitive * P, int i){
	//Get x value, determine U and F based on x. All calculated from initial vars rho, p, v, and gam
	double pos = getx(i);
	if(pos<0.){
		P[i].rho = 1.0;
		P[i].v = 0.0;
		P[i].v_trans = 0.0;
		P[i].W=getW(i,P);
		P[i].p = 1000.0;
		P[i].eps = P[i].p/((gam-1.)*P[i].rho);
		
	}
	else{
		P[i].rho = 1.0;
		P[i].v = 0.0;
		P[i].v_trans = 0.99;
		P[i].W=getW(i,P);
		P[i].p = 0.01;
		P[i].eps = P[i].p/((gam-1.)*P[i].rho);
	}
	U[i].mass=P[i].rho*P[i].W;
	double h = geth(i,P);
	U[i].mom=P[i].rho*h*pow(P[i].W,2)*P[i].v;
	U[i].momT=P[i].rho*h*pow(P[i].W,2)*P[i].v_trans;
	U[i].ene=P[i].rho*h*pow(P[i].W,2)-P[i].p-U[i].mass;

}

void updateF(struct vector * U, struct vector * F, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		F[i].mass = U[i].mass*P[i].v;
		F[i].mom = U[i].mom*P[i].v+P[i].p;
		F[i].momT = U[i].momT*P[i].v;
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
double get_alpha(int i, struct vector * U, struct vector * F, struct primitive * P, int plus){
	double v = .5*(P[i].v+P[i+1].v);
	double cs = .5*(getcs(i,P)+getcs(i+1,P));

	double alpha;
	if(plus==1) alpha = max(0,((v+cs)/(1.+v*cs)));
	else alpha = min(0,((v-cs)/(1.-v*cs)));
	
	//printf("%f", alpha);
	return alpha;	
}

struct vector getFlux(int i, struct vector * U, struct vector * F, struct primitive * P){
	double alpha_plus = get_alpha(i,U,F,P,1);
	double alpha_minus = get_alpha(i,U,F,P,0);
	
	struct vector fhll;
	fhll.mass = (alpha_plus*F[i].mass-alpha_minus*F[i+1].mass+alpha_plus*alpha_minus*(U[i+1].mass-U[i].mass))/(alpha_plus-alpha_minus);
	fhll.mom = (alpha_plus*F[i].mom-alpha_minus*F[i+1].mom+alpha_plus*alpha_minus*(U[i+1].mom-U[i].mom))/(alpha_plus-alpha_minus);
	fhll.momT = (alpha_plus*F[i].momT-alpha_minus*F[i+1].momT+alpha_plus*alpha_minus*(U[i+1].momT-U[i].momT))/(alpha_plus-alpha_minus);
	fhll.ene = (alpha_plus*F[i].ene-alpha_minus*F[i+1].ene+alpha_plus*alpha_minus*(U[i+1].ene-U[i].ene))/(alpha_plus-alpha_minus);
	
	return fhll;
	
}

void updatePrimitive(struct vector * U, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		double delta = 10;
		P[i].p=10.0;
		while(fabs(delta)>.00000001){
			//Sets a minimum value for P so NR doesn't choose an impossible guess
			if(P[i].p<(fabs(U[i].mom)-U[i].ene-U[i].mass)) P[i].p=fabs(U[i].mom)-U[i].ene-U[i].mass;
			P[i].v=U[i].mom/(U[i].ene+U[i].mass+P[i].p);
			P[i].v_trans=U[i].momT/(U[i].ene+U[i].mass+P[i].p);
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
		P[i].v_trans=U[i].momT/(U[i].ene+U[i].mass+P[i].p);
		P[i].W=getW(i,P);
		P[i].rho=U[i].mass/P[i].W;
		P[i].eps=(U[i].ene+U[i].mass*(1.-P[i].W)+P[i].p*(1.+pow(P[i].W,2)))/(U[i].mass*P[i].W);
	
	
	}
}


void advanceTime(struct vector * U, struct vector * F, struct primitive * P){
	struct vector Fiph[N-1];
	int i;
	for(i=0;i<N-1;i++) Fiph[i] = getFlux(i,U,F,P);
	for(i=0;i<N-1;i++){
		U[i].mass -= dt*Fiph[i].mass/dx;
		U[i].mom -= dt*Fiph[i].mom/dx;
		U[i].momT -= dt*Fiph[i].momT/dx;
		U[i].ene -= dt*Fiph[i].ene/dx;
		U[i+1].mass += dt*Fiph[i].mass/dx;
		U[i+1].mom += dt*Fiph[i].mom/dx;
		U[i+1].momT += dt*Fiph[i].momT/dx;
		U[i+1].ene += dt*Fiph[i].ene/dx;
	}
	
	updatePrimitive(U,P); 
	//fixed_bcs(U,F,P);
	initialize(U,P,0);
	initialize(U,P,N-1);
	updateF(U,F,P);

}

		


	
int main(){
	//Loop through different N's, advancing to time t=.25, and output
	N = 500;
	while(N<501){
		dt = time/(8*N);
		dx = length/N;
	
		struct vector U[N];
		struct vector F[N];
		struct primitive P[N];
	
		int i;
		for(i=0;i<N;i++) initialize(U,P,i);
		
		updateF(U,F,P);
	
	
		currentTime=0.;
		while(currentTime<time){
			advanceTime(U,F,P);
			currentTime += dt;
			if(isnan(P[250].rho)==1){
				printf("%f\n",currentTime);
				break;
			}
		}
		
		FILE* f;
		char string[15];
		sprintf(string,"%d",N);
		strcat(string,".txt");
		f = fopen(string,"w");
		
		for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),P[i].rho,P[i].v,P[i].p,P[i].eps);
		//for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),P[i].rho,P[i].W,P[i].p,P[i].v);
		
		
		
		
		fclose(f);
		N = 2*N;
		
	}
	return 0;
}
