/*********************************************
 * OPL 12.10.0.0 Model
 * Author: yue
 * Creation Date: Mar 23, 2022 at 4:42:38 PM
 *********************************************/
int NbParcelles =...; //nombre de parcelles
range Parcelles =1..NbParcelles;
int T =...; //horizon de planification
range Periode=1..T;
int SURF=...;
int lmax= ...; // duree au-delà de laquelle prolonger la jachère n'améliore plus le rendement'
int amax= ...; //nombre max de de semestres de cultures
{int} C[1..2] = ...; //cultures en semestre pair/impair

int s[t in Periode]=(t%2==1)?1:2;

{int} Cultures=C[1] union C[2];

 int Demande[Cultures][Periode]=...;
 
 tuple sommet{
 int l; //age de la jachere
 int a; // age de la culture
 int j; //culture ou jachere en cours 
 }
 
 {sommet} Sommets=...;
 
 tuple arc {
sommet i; //extremite initiale
sommet f; //extremite finale
int rend;
}  
 
 {arc} Arcs=...;
 {arc} InitArc=
 {<<2,0,0>,<2,0,0>,0>,
 <<2,0,0>,<2,1,1>,120>
 };
 
 // variables
 dvar boolean x[Arcs][Periode][Parcelles];
 
/*-----------------------------
* Modelling
*------------------------------*/

// objective
minimize sum(p in Parcelles, a in InitArc) x[a][1][p];
// answer 19, with time 0.37 (s), and 0 explored nodes


subject to {
	// demand constraint
	ctDemand:
		forall(t in Periode, i in C[s[t]])
		  sum(p in Parcelles, a in Arcs : a.f.j==i) x[a,t,p] * a.rend >= Demande[i, t];
 
	// at most one outgoing arc from initial vertex
	ctInitialNode:
		forall( p in Parcelles)
		  sum(a in InitArc) x[a][1][p] <= 1;

	// forbiden arcs
	ctForbiden:
		forall( p in Parcelles)
		  sum(a in Arcs : a not in InitArc) x[a][1][p] == 0;
  
	// conservation flot
	ctConv:
		forall(p in Parcelles, t in 1..T-1, v in Sommets /*: v.j == 0 || v.j in C[s[t]]*/)
		  sum(a in Arcs : a.f==v) x[a][t][p] == sum(a in Arcs : a.i==v) x[a][t+1][p];
}

// Output the current result
execute DISPLAY_RESULT {
  for(c in C[s[2]]){
    writeln(c);
  }
   
}


 