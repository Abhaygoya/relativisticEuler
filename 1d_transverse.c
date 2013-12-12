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
	double ene;
};

struct primitive {
	double rho;
	double v;
	double eps;
	double p;
	double W;
	//Transverse velocity. Does it change?
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

void initialize(struct vector * U, struct primitive * P){
	int i;
	for(i=0; i<N; i++){
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
		double h = 1.+P[i].eps+P[i].p/P[i].rho;
		U[i].mom=P[i].rho*h*pow(P[i].W,2)*P[i].v;
		U[i].ene=P[i].rho*h*pow(P[i].W,2)-P[i].p-U[i].mass;
		
	}
}

void updateF(struct vector * U, struct vector * F, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		F[i].mass = U[i].mass*P[i].v;
		F[i].mom = U[i].mom*P[i].v+P[i].p;
		F[i].ene = U[i].mom-U[i].mass*P[i].v;
	}
}

void fixed_bcs(struct vector * U, struct vector * F, struct primitive * P){
	P[0].rho = 1.0;
	P[0].v = 0.0;
	P[0].W=getW(0,P);
	P[0].p = 1000.0;
	P[0].eps = P[0].p/((gam-1.)*P[0].rho);
	
	U[0].mass=P[0].rho*P[0].W;
	double h = 1.+P[0].eps+P[0].p/P[0].rho;
	U[0].mom=P[0].rho*h*pow(P[0].W,2)*P[0].v;
	U[0].ene=P[0].rho*h*pow(P[0].W,2)-P[0].p-U[0].mass;
	
	P[N-1].rho = 1.0;
	P[N-1].v = 0.0;
	P[N-1].W=getW(N-1,P);
	P[N-1].p = 0.01;
	P[N-1].eps = P[N-1].p/((gam-1.)*P[N-1].rho);
	
	U[N-1].mass=P[N-1].rho*P[N-1].W;
	h = 1.+P[N-1].eps+P[N-1].p/P[N-1].rho;
	U[N-1].mom=P[N-1].rho*h*pow(P[N-1].W,2)*P[N-1].v;
	U[N-1].ene=P[N-1].rho*h*pow(P[N-1].W,2)-P[N-1].p-U[N-1].mass;
	
	
	
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
	fhll.ene = (alpha_plus*F[i].ene-alpha_minus*F[i+1].ene+alpha_plus*alpha_minus*(U[i+1].ene-U[i].ene))/(alpha_plus-alpha_minus);
	
	return fhll;
	
}

void updatePrimitive(struct vector * U, struct primitive * P){
	int i;
	for(i=0;i<N;i++){
		double delta = 10;
		P[i].p=10.0;
		while(fabs(delta)>.0000000000001){
			//Sets a minimum value for P so NR doesn't choose an impossible guess
			if(P[i].p<(fabs(U[i].mom)-U[i].ene-U[i].mass)) P[i].p=fabs(U[i].mom)-U[i].ene-U[i].mass;
			P[i].v=U[i].mom/(U[i].ene+U[i].mass+P[i].p);
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
		U[i].ene -= dt*Fiph[i].ene/dx;
		U[i+1].mass += dt*Fiph[i].mass/dx;
		U[i+1].mom += dt*Fiph[i].mom/dx;
		U[i+1].ene += dt*Fiph[i].ene/dx;
	}
	
	updatePrimitive(U,P); 
	fixed_bcs(U,F,P);
	updateF(U,F,P);

}

		


	
int main(){
	//Loop through different N's, advancing to time t=.25, and output
	N = 500;
	while(N<501){
		dt = time/(10*N);
		dx = length/N;
	
		struct vector U[N];
		struct vector F[N];
		struct primitive P[N];
	
		initialize(U,P);
		updateF(U,F,P);
	
	
		currentTime=0.;
		while(currentTime<.000001){
			advanceTime(U,F,P);
			currentTime += dt;
		}
		
		FILE* f;
		char string[15];
		sprintf(string,"%d",N);
		strcat(string,".txt");
		f = fopen(string,"w");
		
		int i;
		for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),P[i].rho,P[i].v,P[i].p,P[i].eps);
		//for(i=0;i<N;i++) fprintf(f,"%15.15f     %15.15f       %15.15f    %15.15f    %15.15f\n", getx(i),P[i].rho,P[i].W,P[i].p,P[i].v);
		
		
		
		
		fclose(f);
		N = 2*N;
		
	}
	return 0;
}
