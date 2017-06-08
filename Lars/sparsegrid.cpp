#include <stdio.h>
#include <stdlib.h>
#include <math.h>




/* Folgendes ist die Ausgangssituation: Algorithmus 3 (enumeration) liefert Vektoren 
(oder wie man es nennen will), deren Eintraege die Levels von den einzelnen
 eindimensionalen Quadraturregeln sind. 

Von diesen verschiedenen eindimensionalen Levels müssen dann alle moeglichen 
Kombinationen gebildet werden. (aehnlich wie beim product grid) 
Der Unterschied zum Product grid ist, dass die Vektoren meistens verschiedene Laengen haben.
Zum Bilden all dieser Kombinationen dient die Funktion enumeration2.
*/

// Funktion für nodes von clenshaw curtis
double* clen_curtn(int l)
{
	int Nl = pow(2,l)-1;
	double* nodes = new double[Nl];
	for(int i=1; i<=Nl; i++)
	{
		nodes[i-1]=(.5*(1.-cos((double)(M_PI*i/(Nl+1.)))));
	}
	return nodes;
}


int enumeration(int *k, int d, int l){
	int S_k=d;
	while(1<2){

// Bis hier ist es der Algorithmus abgeschrieben vom Arbeitsblatt, ab hier kommt etwas anderes
// 

	// Das hier dient nur zur Demonstration, dass es überhaupt klappt. So habe ich heute Mittag
	// auch die Bilder fuer Clenshaw Curtis erstellt.
	
	// k[0] und k[1] enthalten die Levels im zweidimensionalen sparse grid
	double* nodes1 = clen_curtn(k[0]);
	double* nodes2 = clen_curtn(k[1]);
	
	// Hier werden einfach nur die Stuetzstellen in eine Datei geschrieben
	// Der Algorithmus passiert diese Stelle mehrmals, und jedes mal wird die Datei geoeffnet
	// und die zusaetzlichen Stuetzstellen werden eingetragen. Mit gnuplot kann man daraus dann 
	// diese Gitter erzeugen.
	
	FILE *fp;
	fp = fopen("combinations", "a");
	
	for(int i = 0; i<pow(2,k[0])-1; i++){
		for(int j=0; j<pow(2,k[1])-1; j++){
			fprintf(fp,"%f  %f\n", nodes1[i], nodes2[j]);
		}
	}
	fclose(fp);
	
	


	// zur Kontrolle, damit man sehen kann, welche Kombinationen der Algorithmus erzeugt
	for(int i=0; i<d; i++){
		printf("%i \n", k[i]);
	}
	printf("end \n"); 
	
	// Hier ist es wieder einfach der Algorithmus vom Arbeitsblatt
	
		for(int j=1; j<=d; j++){
			k[j-1]=k[j-1]+1;
			S_k = S_k +1;
			if(S_k > d+l-1){
				if(j==d){ return 0; }
				S_k = S_k-k[j-1]+1;
				k[j-1] = 1;
			}
			else break;
		}
		
	}
	
	return 1;
}

/* Die obige Funktion hat das Problem, dass sie wegen der verschachtelten For-Schleife
nur fuer zweidimensionale sparse grids funktioniert. Der Algorithmus enumeration kann allerdings
mit beliebigen Dimensionen umgehen.
enumeration erzeugt also im zweidimensionalen Fall Paare von Zahlen, die angeben welche Paare
von Levels zu dem Sparse Grid dazugehören sollen. Das sind alle Grids, die auf Seite 7 in Fig.1
oberhalb der durchgezogenen Trennlinie liegen (und fuer die k_1+...+k_d<=l+d-1 gilt).

Im allgemeinen Fall ist jetzt das Problem, dass man aus diesen Vektoren verschiedener Laenge
Punkte bilden muss. Diese Punkte kann man dann nutzen, um sie in der mehrdimensionalen Funktion
einzusetzen oder um das Produkt der eindimensionalen Gewichte zu bestimmen. 
Diese Punkte werden vom folgenden Algorithmus erzeugt. */

void enumeration2(int* vec, int* klevel, int d){
	int I;
	int count;
	int last = d-1;
	
	while(last+1==d){
	
	// Hier werden die Kombinationen, die der Algorithmus gerade erzeugt, ausgegeben
	for(int j=0; j<d; j++){
		printf("  %i", vec[j]);
	}
	printf(" end \n");
	
		if(vec[last]<klevel[last]){
			vec[last]++;
			}
		else {
		for(int i=d-1; i>=0; i--){
			if(vec[i]<klevel[i]) {
			vec[i]++;
			break;
			}
			
			count = 0;
			for(I=0; I<d; I++){
			if(vec[I]==klevel[I]) { count++; 
				}
			}
			if(count==d) d=-2; 
			
			if(vec[i]==klevel[i]) vec[i]=1;
			}

		}
	}
}

/* enumeration2 braucht als Input einen Vektor (1,1,..,1), einen Array klevel, der die Anzahl 
der Stutzstellen/Gewichte in dem Level enthaelt. Also ist klevel[i]=pow(2,k[i])-1, mit k[i]
wie auf dem Arbeitsblatt. d ist die Dimension und damit die Laenge von klevel.
Der Algorithmus erzeugt die Vektoren in folgender Weise:

klevel=[3,3,3], das heisst, k[i]=2 und auf jedem Level gibt es drei Stuetzstellen/Gewichte
d = 3

1 1 1
1 1 2
1 1 3
1 2 1
1 2 2 
1 2 3
1 3 1
1 3 2 
1 3 3
2 1 1
2 1 2
usw...

bis man in diesem Falle 3*3*3=27 Kombinationen erreicht hat und am Ende 3 3 3 steht.
Ich bin selber etwas erstaunt, dass es so geklappt hat, aber ich habe den Algorithmus fuer 
viele verschiedene Sachen getestet und er hat sich als voll funktionsfaehig erwiesen.
Natuerlich kann klevel auch anders aussehen, z.b. klevel=[3,3,7,1,7,3]

Genau so wie beim enumeration Algorithmus vom Arbeitsblatt werden im Algorithmus alle Kombinationen 
erzeugt, zurueckgegeben wird eigentlich nichts.
*/

/* Was man im ersten Teil sehen konnte war, wie Nodes erstellt werden. An der Stelle haette 
man mit enumeration2 auch die Verwendung der verschachtelten For-Schleife umgehen koennen.

Der Folgende Teil war ein Versuch von mir, moeglichst allgemein die Berechnung der Gewichte 
durchzufuehren, die ja die Multiplikation der eindimensionalen Gewichte aus allen Levels
erfordert. Das ist noch sehr unuebersichtlich und erforderte die Veraenderung der Funktion 
enumeration2, da inmitten der Funktion nun immer andere Funktionen aufgerufen werden. 
*/

void trap_rulew(int l, double* weights){
	int Nl = pow(2, l)-1;
	
	for(int i=1; i<=Nl; i++){
		 weights[i-1] = (double) 1/(Nl+1);
		
	}
	weights[0] = (double) 3/(2*(Nl+1));
	weights[Nl-1] = (double) 3/(2*(Nl+1));
}


void enumeration_weights(int* vec, int* klevel, int d, int last, double* weights, int* k){
	int I;
	int count;
	int c2=0;
	int l=0;
	double product = 1;
	
	// sum_ki soll die Laenge des Vektors mit allen Gewichten aus allen Levels sein
	int sum_ki = 10;
	double* ccweights = new double[sum_ki];
	double* cctemp = new double[sum_ki];
	
	// Umstaendlicher shit auf den ich an dieser Stelle nicht weiter eingehen moechte, 
	// machen wir noch anders, klappt aber so schon so ungefaehr
	trap_rulew(k[0], ccweights);	
	for(int i=1; i<d; i++){
		trap_rulew(k[i], cctemp);
		for(int jj=0; jj<klevel[i]; jj++){
			ccweights[jj+klevel[i-1]]=cctemp[jj];
		}
	}
	
	// Hier werden einmal die Gewichte vor dem Bilden der d-dimensionalen Produkte ausgegeben
	for(int i=0; i<sum_ki; i++){
		printf("%f weights \n", ccweights[i] );
	}
	
	while(last+1==d){
	
	
	// hier wird dann das d-dimensionale Produkt gebildet und im Terminal ausgegeben
	l=0;
	for(int j=0; j<d; j++){
		klevel[-1]=0;
		product = ccweights[klevel[j-1]+vec[j]-1]*product;
		l++;
	}
	
	weights[c2] = product;
	c2++;
	product = 1;
	
	
	
		if(vec[last]<klevel[last]){
			vec[last]++;
			}
		else {
		for(int i=d-1; i>=0; i--){
			if(vec[i]<klevel[i]) {
			vec[i]++;
			break;
			}
			
			count = 0;
			for(I=0; I<d; I++){
			if(vec[I]==klevel[I]) { count++; 
				}
			}
			if(count==d) d=-2; 
			
			if(vec[i]==klevel[i]) vec[i]=1;
			}

		}
	}
}







int main () {
	
	// l ist das Level, d die Dimension, alles wie auf dem Arbeitsblatt
	int l = 3;
	int d = 2;
	int* k = new int[d]; 
	
	// Das ist dafür da, dass enumeration mit k=1 anfangen kann
	for(int i=0; i<d; i++){
		k[i] = 1;
	}
	
	enumeration(k, d, l);
	
	// bis hierhin war alles fuer enumeration, ab hier ist alles fuer enumeration_weights
	
	int* klevel = new int[6];
	int* vec = new int[6];
	double* weights = new double[21];

	for(int i=0; i<6; i++){
		vec[i]=1;
	}
	
	klevel[0]=3;
	klevel[1]=3;
	klevel[2]=3;
	
	enumeration_weights(vec, klevel, 3, 2, weights, k);

	// Die Ausgabe der Gewichte, die Summe sollte 1 ergeben (und das tut sie auch)
	for(int i=0; i<27; i++){
	printf("%f+", weights[i]);
	}
	
	/* Zusammenfassung:
	Was passiert ist Folgendes: Man entscheidet sich fuer eine Dimension und ein Level.
	Auf Seite 7 ist beschrieben, wie dann daraus ein sparse grid erstellt wird.
	Das geschieht im Algorithmus durch enumeration. 
	Wenn man so ein sparse grid hat, muss man aber noch die Gewichte berechnen.
	Dafuer kann man enumeration nicht nutzen, denn enumeration hat nur angegeben, 
	welche Paare von Levels fuer das sparse grid zulaessig sind. 
	enumeration2 nimmt nochmal jedes einzelne dieser Paare und bildet alle Kombinationen.
	Wenn man aus enumeration erhaelt, dass (2,3) ein zulaessiges Paar fuer das Sparse grid ist,
	dann uebergibt man diese Information an enumeration2. enumeration2 kann dann auf dem 
	Gitter aus den zwei Vektoren der Laenge 2^2-1=3 und 2^3-1=7 alle Punkte dieses
	Gitters bestimmen. Diese Punkte (i,j) sind dann das i-te Gewicht aus k_1 und das 
	j-te Gewicht aus k_2. So kann man dann diese Werte Multiplizieren, sodass man auf die 
	zweidimensionalen Gewichte kommt. Das fuehrt dazu, dass man auf dem sparse grid alle
	Stuetzstellen und alle dazugehoerigen Gewichte erhaelt, und dann setzt man die Stuetzstellen
	in f ein, multipliziert an allen Stellen mit den Gewichten und summiert dann ueber alles
	und man hat den Wert des Integrals. Alles ganz einfach also.
	*/
	
}



