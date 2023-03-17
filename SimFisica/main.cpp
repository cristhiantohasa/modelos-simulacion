#include "cabecera.h"
#include "MathFunctions.h"

using namespace std;

//---------------------------------------------------------------------------
room        r;                                  //Cubic Room
source      s1;
receptor    r1;
double      tiempo,distancia;
bool        cuartoCargado,TestExtra,Recorrido;
reflection  *ry;
int         numRay,eNumRay;                     //Numero de rayos que usar� la simulacion

//Visualizacion del ray tracing
vector      **vec_gra;          //Vectores para grafico (por rayo y por reflexi�n)
double      **mod_vg;           //M�dulo acumulativo (por rayo y por reflexi�n)
double      **tim_vg;           //Tiempo transcurrido entre cada reflexi�n
double      v_son;              //Velocidad del sonido (ej: 340 m/s)
double      eneRay;             //Energ�a de un rayo
double      alfa, delta;        //alfa es coeficiente de absorci�n y delta es coeficiente de difusi�n
matEnergia  mE,mEr;             //Matriz de energ�a distribuida en el espacio y tiempo
matDouble   mD,mDr, mP, mPr;    //Matriz double y de porcentajes

//Punto 3
matTiempo   mT,mTr;     //Matriz double
int         NTRI;       //N�mero total de tri�ngulos que tiene la sala
//---------------------------------------------------------------------------

void CargarCubo();
void Simulacion(int x0, int y0, int z0);

int main() {

    cuartoCargado=false;
    TestExtra=false;
    Recorrido=false;
    numRay=0;
    ry=NULL;
    eNumRay=1000;

    int caso=-1;
    bool nosalir=true;

    while(nosalir){

        cout<<"MENU\n";
        cout<<"1. Cargar sala y simular.\n";
        cout<<"2. Salir.\n";
        cout<<"Escoja una opcion:\t";
        cin>>caso;

        switch(caso){
            case 1:
                CargarCubo();
                cout<<"\nSala cargada.\n\n";
                for(int i = -1 ; i <= 1 ; i++)
                    for (int j = -1 ; j <= 1 ; j++)
                        for ( int k = -1 ; k <= 1 ; k++ ){
                            TestExtra=false;
                            Simulacion(i,j,k);
                        }
                cout<<"\nSimulacion realizada.\n\n";
                break;
            case 2:
                nosalir=false;
                cout<<"Adios.\n\n";
                break;
            default:
                cout<<"Opcion invalida.\n\n";
                break;
        }

    }

    return 0;

}

//---------------------------------------------------------------------------
void CargarCubo(){

    int nt;     //Contador de tri�ngulos para el ID �nico
    int nd=2;   //Numero de divisiones de las caras cuadradas

    cuartoCargado=true;

    //Inicializaciones cubo
    r.Clear();
    r.NewPlanes(6);
    //-------------square 1 back
    r.p[0].NewPoints(4);
    r.p[0].p[0].x=-2.0f;
    r.p[0].p[0].y=2.0f;
    r.p[0].p[0].z=2.0f;
    r.p[0].p[1].x=-2.0f;
    r.p[0].p[1].y=-2.0f;
    r.p[0].p[1].z=2.0f;
    r.p[0].p[2].x=-2.0f;
    r.p[0].p[2].y=-2.0f;
    r.p[0].p[2].z=-2.0f;
    r.p[0].p[3].x=-2.0f;
    r.p[0].p[3].y=2.0f;
    r.p[0].p[3].z=-2.0f;
    r.p[0].MoreTriangle(nd);    //PointGenTriangle();

    //-------------square 2 front
    r.p[1].NewPoints(4);
    r.p[1].p[0].x=2.0f;
    r.p[1].p[0].y=2.0f;
    r.p[1].p[0].z=2.0f;
    r.p[1].p[1].x=2.0f;
    r.p[1].p[1].y=2.0f;
    r.p[1].p[1].z=-2.0f;
    r.p[1].p[2].x=2.0f;
    r.p[1].p[2].y=-2.0f;
    r.p[1].p[2].z=-2.0f;
    r.p[1].p[3].x=2.0f;
    r.p[1].p[3].y=-2.0f;
    r.p[1].p[3].z=2.0f;
    r.p[1].MoreTriangle(nd);    //PointGenTriangle();

    //-------------square 3 left
    r.p[2].NewPoints(4);
    r.p[2].p[0].x=-2.0f;
    r.p[2].p[0].y=-2.0f;
    r.p[2].p[0].z=2.0f;
    r.p[2].p[1].x=2.0f;
    r.p[2].p[1].y=-2.0f;
    r.p[2].p[1].z=2.0f;
    r.p[2].p[2].x=2.0f;
    r.p[2].p[2].y=-2.0f;
    r.p[2].p[2].z=-2.0f;
    r.p[2].p[3].x=-2.0f;
    r.p[2].p[3].y=-2.0f;
    r.p[2].p[3].z=-2.0f;
    r.p[2].MoreTriangle(nd);    //PointGenTriangle();

    //-------------square 4 right
    r.p[3].NewPoints(4);
    r.p[3].p[0].x=2.0f;
    r.p[3].p[0].y=2.0f;
    r.p[3].p[0].z=2.0f;
    r.p[3].p[1].x=-2.0f;
    r.p[3].p[1].y=2.0f;
    r.p[3].p[1].z=2.0f;
    r.p[3].p[2].x=-2.0f;
    r.p[3].p[2].y=2.0f;
    r.p[3].p[2].z=-2.0f;
    r.p[3].p[3].x=2.0f;
    r.p[3].p[3].y=2.0f;
    r.p[3].p[3].z=-2.0f;
    r.p[3].MoreTriangle(nd);    //PointGenTriangle();

    //-------------square 5 top
    r.p[4].NewPoints(4);
    r.p[4].p[0].x=-2.0f;
    r.p[4].p[0].y=-2.0f;
    r.p[4].p[0].z=2.0f;
    r.p[4].p[1].x=-2.0f;
    r.p[4].p[1].y=2.0f;
    r.p[4].p[1].z=2.0f;
    r.p[4].p[2].x=2.0f;
    r.p[4].p[2].y=2.0f;
    r.p[4].p[2].z=2.0f;
    r.p[4].p[3].x=2.0f;
    r.p[4].p[3].y=-2.0f;
    r.p[4].p[3].z=2.0f;
    r.p[4].MoreTriangle(nd);    //PointGenTriangle();

    //-------------square 1 bottom
    r.p[5].NewPoints(4);
    r.p[5].p[0].x=-2.0f;
    r.p[5].p[0].y=2.0f;
    r.p[5].p[0].z=-2.0f;
    r.p[5].p[1].x=-2.0f;
    r.p[5].p[1].y=-2.0f;
    r.p[5].p[1].z=-2.0f;
    r.p[5].p[2].x=2.0f;
    r.p[5].p[2].y=-2.0f;
    r.p[5].p[2].z=-2.0f;
    r.p[5].p[3].x=2.0f;
    r.p[5].p[3].y=2.0f;
    r.p[5].p[3].z=-2.0f;
    r.p[5].MoreTriangle(nd);    //PointGenTriangle();

    //Inicializar normales de los planos
    //Se calculan las normales con la normal de su primer triangulo
    nt=0;
    for (int i=0; i<r.NP; i++){

        r.p[i].n=TriangleNormal(r.p[i].t[0]);

        for (int j=0; j<r.p[i].NT; j++){

            r.p[i].t[j].CalculateProjection();
            r.p[i].t[j].calcularBC();           //baricentro
            r.p[i].t[j].ID=nt;                  //identificador �nico de cada tri�ngulo de la sala

            nt++;

        }

    }

    NTRI=nt;
    r.MaxDistance();
    mD.init(NTRI, NTRI);
    mDr.init(NTRI,1);

    //Punto 3
    mT.init(NTRI, NTRI);
    mTr.init(NTRI,1);

    //punto 4
    mP.init(NTRI, NTRI);
    mPr.init(NTRI,1);

    //punto 5
    double * areaTotal= new double[NTRI];
    for (int i=0; i<NTRI; i++){
        areaTotal[i]=0.0;
    }

    for (int i=0; i<r.NP; i++){                     //planos
        for (int j=0; j<r.p[i].NT; j++){            //triangulos
            for (int k=0; k<r.NP; k++){             //planos
                for (int l=0; l<r.p[k].NT; l++){    //triangulos
                    if(i!=k){

                        mD.d[r.p[i].t[j].ID][r.p[k].t[l].ID]=r.p[i].t[j].bc.distancia(r.p[k].t[l].bc);      //distancia
                        mDr.d[r.p[i].t[j].ID][0]=r.p[i].t[j].bc.distancia(r1.p);                            //distancia receptor

                        //Punto 3
                        mT.t[r.p[i].t[j].ID][r.p[k].t[l].ID]= mD.d[r.p[i].t[j].ID][r.p[k].t[l].ID]/V_SON;   //tiempo
                        mTr.t[r.p[i].t[j].ID][0]= mD.d[r.p[i].t[j].ID][0]/V_SON;                            //tiempo  receptor

                        //Punto 4
                        mP.d[r.p[i].t[j].ID][r.p[k].t[l].ID]=r.p[i].t[j].calcularAngSolido(r.p[k].t[l]);    //areas
                        mPr.d[r.p[i].t[j].ID][0]=r1.calcularAnguloSolido(r.p[i].t[j].bc);                   //areas
                        areaTotal[r.p[i].t[j].ID]+=mP.d[r.p[i].t[j].ID][r.p[k].t[l].ID];                    //suma las �reas peque�as

                    }
                }
            }
        }
    }

    for(int i=0; i<NTRI; i++){
        for(int j=0; j<NTRI; j++){
            mP.d[i][j]=mP.d[i][j]/areaTotal[i]; //porcentaje
        }
        mPr.d[i][0]=mPr.d[i][0]/areaTotal[i];
    }

    mD.guardarArchivo(nt, "Distancia");
    mDr.guardarArchivoR(nt, "DistanceReceptor");
    mT.guardarArchivo(nt);
    mTr.guardarArchivoR(nt);
    mP.guardarArchivo(nt, "Angulo");
    mPr.guardarArchivoR(nt,"AnguloReceptor");

}
//---------------------------------------------------------------------------
void Simulacion(int x0, int y0, int z0){

    if(TestExtra==false){

        if(s1.NRAYS!=0){
            delete ry;
            ry=NULL;
        }

        //Definimos la duraci�n de la simulacion en milisegundos
        int durSim=1000;

        //Definimos la posici�n de la fuente s1
        s1.p.x=1.5;
        s1.p.y=1.5;
        s1.p.z=1.5;

        //Definimos la posici�n del receptor r1
        r1.p.x=x0;   //-1 0 1
        r1.p.y=y0;   //-1 0 1
        r1.p.z=z0;   //-1 0 1

        //Definimos energia de la fuente
        s1.eF=100;

        //Definimos los espacios de tiempo en el receptor
        r1.createTimeSamples(durSim);

        //Definimos coeficientes de absorci�n y difusi�n
        alfa=0.2;
        delta=0.15;

        //Velocidad del sonido
        v_son=340;

        //Tomo el numero de rayos que solicite el usuario
        numRay=eNumRay;

        //Defino la direccion de los rayos para la fuente s1
        s1.createRays(numRay);

        //Ajusto la cantidad de rayos a lo que determina el procedimiento de subdivisi�n del icosaedro
        numRay=s1.NRAYS;
        eNumRay=numRay;

        //Me aseguro que no haya nada en ry
        delete ry;

        ry=r.RayTracing(s1.p,s1.Rays,s1.NRAYS);

        eneRay=s1.eF/s1.NRAYS;

        if (mE.tri!=0 || mE.tim!=0){
            mE.clear();
        }

        mE.init(NTRI,durSim);

        //Variable para la visualizacion grafica
        vec_gra = new vector*[numRay];
        mod_vg = new double*[numRay];
        tim_vg = new double*[numRay];

        for(int R=0; R<numRay; R++){    //Recorro todos los rayos que parten de la fuente

            double eneRes=eneRay;
            vec_gra[R] = new vector[ry[R].N-1];
            mod_vg[R] = new double[ry[R].N];
            mod_vg[R][0]=0;
            tim_vg[R] = new double[ry[R].N];
            tim_vg[R][0]=0;

            for (int i=0; i<ry[R].N-1; i++){    //genera el vector a ser reflejado y luego calcula su m�dulo
                vec_gra[R][i]=ry[R].r[i+1]-ry[R].r[i];
                mod_vg[R][i+1]=Module(vec_gra[R][i]);
                tim_vg[R][i+1]=1000*(mod_vg[R][i+1]/v_son); //Calcula el tiempo en ms
            }

            for (int i=1; i<ry[R].N; i++){

                mod_vg[R][i]=mod_vg[R][i-1]+mod_vg[R][i];
                tim_vg[R][i]=tim_vg[R][i-1]+tim_vg[R][i];

                int indTim, indTri;

                indTim = int(tim_vg[R][i]);
                indTri = ry[R].idTriangle[i-1];

                if(i<ry[R].N-1){

                    double dis;
                    int tim;
                    vector vec_nr,vec_uni;
                    point p_r,p_o;
                    vec_uni=Normal(vec_gra[R][i]);
                    vec_nr=vec_uni*(-1);
                    p_r=r1.p;
                    p_o=ry[R].r[i-1];
                    dis=IntersectionDistance(vec_nr,p_r,vec_uni,p_o);

                    if(dis>0&&dis<r.maxd){
                        p_r=p_o+(vec_uni*dis);
                        if(Module(p_r-r1.p)<r1.ReceptionRadius){
                            dis=Module(vec_uni*dis);
                            tim=int(1000*dis/v_son);
                            r1.eR[indTim+tim]+=eneRes;
                        }
                    }
                }

                mE.ene[indTri][indTim]+=eneRes*(1-alfa)*delta; //energia difusa del rayo que se almacena en la matriz de energia;
                eneRes=eneRes*(1-alfa)*(1-delta);

            }

        }

        //Paso 5
        for(int i=0; i<durSim; i++){
            for(int j=0; j<NTRI; j++){  //triangulo de partida
                if(mE.ene[j][i]!=0){
                    for(int k=0; k<NTRI; k++){  //tri�ngulo de destino
                        if(i+mT.t[j][k]<durSim){
                            mE.ene[k][i+mT.t[j][k]]+=mE.ene[j][i]*mP.d[j][k]*(1-alfa);
                        }
                    }
                }
                if(i+mTr.t[j][0]<durSim){
                    r1.eR[i+mTr.t[j][0]]+=mE.ene[j][i]*mPr.d[j][0];
                }
            }
        }

        mE.guardarArchivo("Energia");

        //paso 6
        r1.guardarArchivo();
        TestExtra=true;

    }
    else
    {
        TestExtra=false;
    }
}
//---------------------------------------------------------------------------
