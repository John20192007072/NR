clc;clear;
%sad
%e=input('ingrese el error maximo');
%iterMAX=input('Ingrese el número de iteraciones maximas para lograr el error');
% nodos=input('Ingrese número de nodos:');
% lineas=input('Ingrese número de ramas:');
% L=zeros(lineas,5);
% N=zeros(nodos,8);
% for K=1:nodos
%     N(K,1)=K;
% end
% for K=1:nodos
%     fprintf(['ingrese que tipo de nodo (slack=0 PV=1 PQ=2) es el nodo %i\n'],K);
%     N(K,2)=input('');
% end
% for K=1:nodos
% if N(K,2)==0
%     fprintf('ingrese la magnitud de la tension en el nodo %i:',K);
%     N(K,7)=input('');
%     fprintf('ingrese el angulo de la tension en el nodo %i:',K);
%     N(K,8)=input('')
% end 
% if N(K,2)==1
%     fprintf('ingrese la magnitud de la tension en el nodo %i:',K);
%     N(K,7)=input('');
%     fprintf('ingrese la potencia activa demandada en el nodo %i:',K);
%     N(K,5)=input('');
%     fprintf('ingrese la potencia reactiva demandada en el nodo %i:',K);
%     N(K,6)=input('');
%     fprintf('ingrese la potencia activa generada en el nodo %i:',K);
%     N(K,3)=input('');
%     fprintf('ingrese la potencia reactiva generada en el nodo%i:',K);
%     N(K,4)=input('');
% end 
% if N(K,2)==2
%     fprintf('ingrese la potencia activa demandada en el nodo %i:',K);
%     N(K,5)=input('');
%     fprintf('ingrese la potencia reactiva demandada en el nodo %i:',K);
%     N(K,6)=input('');
%     fprintf('ingrese la potencia activa generada (generacion distribuida) en el nodo %i:',K);
%     N(K,3)=input('');
%     fprintf('ingrese la potencia reactiva generada (generacion distribuida) en el nodo%i:',K);
%     N(K,4)=input('');
% end 
%     fprintf('ingrese la potencia activa generada en el nodo %i:',K);
%     N(K,3)=input('');
%     fprintf('ingrese la potencia reactiva generada en el nodo%i:',K);
%     N(K,4)=input('');
%     fprintf('ingrese la potencia activa demandada en el nodo %i:',K);
%     N(K,5)=input('');
%     fprintf('ingrese la potencia reactiva demandada en el nodo %i:',K);
%     N(K,6)=input('');
%     fprintf('ingrese la magnitud de la tension en el nodo %i:',K);
%     N(K,7)=input('');
%     fprintf('ingrese el angulo de la tension en el nodo %i:',K);
%     N(K,8)=input('');  
% end
% 
% for k=1:lineas
% fprintf('para la rama %i ingrese los nodos donde se conecta \n',k);
% fprintf('ingrese el nodo inicial:');
% L(k,1)=input('');
% fprintf('ingrese el nodo Final:');
% L(k,2)=input('');
% end
% for k=1:lineas
% fprintf('para la rama %i ingrese la conductancía\n',k);
% L(k,3)=input('');
% fprintf('para la rama %i ingrese la suceptancia\n',k);
% L(k,4)=input('');
% fprintf('para la rama %i ingrese el valor de Yc/2 \n',k);
% L(k,5)=input('');
% end
% L
Sbase=100;%MVA
VB=230;%KV
%[K M   Rkm     Xkm     Yc/2]
L=[1 2 0.02 0.04 0;
    1 3 0.01 0.03 0;
    2 3 0.0125 0.025 0;];
%T slack=0 PV=1 PQ=2
%[K T Pgk Qgk Pck Qck    V    A]
N=[1 0 0   0   0 0 1.05    0;
    2 2 0   0   4 2.5 1    0;
    3 1 2   0   0 0 1.04    0;];
%pgloo crow

FL= size(L,1);
FN= size(N,1);
Ybarra=zeros(FN,FN);
for j=1:FL
    K=L(j,1);M=L(j,2);Zkm=L(j,3)+1j*L(j,4);Yc=1j*L(j,5);
    Ybarra(K,M)=-1/Zkm;%Triangular Superior
    Ybarra(M,K)=-1/Zkm;%Triangular Inferior
    Ybarra(K,K)=Ybarra(K,K)+1/Zkm+Yc;%Elementos dte la diagonal
    Ybarra(M,M)=Ybarra(M,M)+1/Zkm+Yc;%Elementos de la diagonal
end
U=N
wy=N;
V0=N(:,7);
VA=N(:,8);
VAN=VA;
V0N=V0;

Ns=0;
Qs=0;
error= input('Ingrese el error en los calculos ');
iterMAX= input('Ingrese el número de iteraciones maximas para lograr el error: ');
for X=1:iterMAX
%Obteniendo la potencia activa y reactiva
FP=0;        % # Filas de la matriz de potencias calculadas
NPQ=0;       % # De nodos PQ
NPV=0;       % # De nodos PV
for j=1 :FN  % Se crea la matriz de potencias y su dimencion depende de el # de nodos pv y pq
    H=N(j,2);FP=H+FP;
    if H==2
        NPQ=NPQ+1;
    end
    if H==1
        NPV=NPV+1;
    end
end
magnitud=(NPQ+NPV);
P=zeros(FP,1);
Pn=zeros(FP,1);
Po=zeros(FP,1);
JA=zeros(FP,FP);
Vant=zeros(FP,1);%% usada para capturar los valores angulo y tencion en los nodos PV y PQ de la matriz N

Pcalculadas=zeros(magnitud,1);
Qcalculadas=zeros(magnitud,1);
%Creando matrizes para almacenar los valores actualizados de magnitud y
%angulo de tension
%CREANDO VECTOR DE POTENCIAS INICIALES
V=1; 
for j=1:FN
    if N(j,2)==2||N(j,2)==1
        Po(V,1)=N(j,3)-N(j,5);
        V=V+1;
        
    end
end
for j=1:FN
    if N(j,2)==2
        Po(V,1)=N(j,4)-N(j,6);
        V=V+1;
    end
end
%_______________________________
i=1;% Contador para llevar la cuenta de los numeros que colocan en la matriz P
F=1; %Controla la posicion de las filas del jacobiano
m=1;

for j=1:FN

        if N(j,2)==1 || N(j,2)==2 
       for t=1:FN  
  P(m,1)=(-1)*(N(j,7)*N(t,7)*abs(Ybarra(j,t))*sin(N(t,8)-N(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qcalculadas(m,1)=P(m,1)+Qcalculadas(m,1);
       end
       
       m=m+1;
        end
         

      if N(j,2)==1 || N(j,2)==2 
  for t=1:FN
       P(i,1)=N(j,7)*N(t,7)*abs(Ybarra(j,t))*cos(N(t,8)-N(j,8)+angle(Ybarra(j,t)));
  Pcalculadas(i,1)=P(i,1)+Pcalculadas(i,1);
  end
      end
      
    if N(j,2)==2||N(j,2)==1
        for t=1:FN
        P(i,1)=N(j,7)*N(t,7)*abs(Ybarra(j,t))*cos(N(t,8)-N(j,8)+angle(Ybarra(j,t)));% hallando P calculada. El valor de p calculada queda en Pn y no en P
        Pn(i,1)=P(i,1)+Pn(i,1);
        end
        i=i+1;
    end
        Q=1;   %Controla la posicion de las columnas del jacobiano
        if N(j,2)==2||N(j,2)==1
        for t=1:FN   
            if N(t,2)==2||N(t,2)==1
                if j==t
                    JA(F,Q)=-Qcalculadas(F,1)-(imag(Ybarra(j,j))*(N(j,7))^2); %Calcula valores de la diagonal de la matriz H 
                end
                if j~=t
                   JA(F,Q)=N(j,7)*N(t,7)*(real(Ybarra(j,t))*sin(N(j,8)-N(t,8))-imag(Ybarra(j,t))*cos(N(j,8)-N(t,8))); %Calcula valores de la triangula superior y inferior de la matriz H 
                end                
            Q=1+Q;
             end
            end
            for t=1:FN 
                if N(t,2)==2
                if j==t

                    JA(F,Q)=Pcalculadas(F,1)+real(Ybarra(j,j))*(N(j,7)^2); %Calcula valores de la diagonal de la matriz N crow
                end
                if j~=t
                    JA(F,Q)=N(j,7)*N(t,7)*((real(Ybarra(j,t)))*cos(N(j,8)-N(t,8))+imag(Ybarra(j,t)*sin(N(j,8)-N(t,8)))); %Calcula valores de la triangula superior y inferior de la matriz N crow

                end 
                    Q=1+Q;
             end
            end
           F=F+1;
        end 
        
end
%_______________________________

for j=1:FN
        Q=1;   %Controla la posicion de las columnas del jacobiano
        
        if N(j,2)==2
            o=1;
        for t=1:FN   
            if N(t,2)==2||N(t,2)==1

                
                if j==t
                    JA(F,Q)=Pcalculadas(o,1)-real(Ybarra(j,j))*(N(j,7)^2); %Calcula valores de la diagonal de la matriz J crow
                end
                  o=o+1;
                if j~=t
                    JA(F,Q)=-N(j,7)*N(t,7)*((real(Ybarra(j,t)))*cos(N(j,8)-N(t,8))+imag(Ybarra(j,t)*sin(N(j,8)-N(t,8)))); %Calcula valores de la triangula superior y inferior de la matriz J crow
                end
            Q=1+Q;
         
            end
            
        end

        lt=1;

            for t=1:FN 
                if N(t,2)==2
                if j==t
                    JA(F,Q)=Qcalculadas(lt,1)-(imag(Ybarra(j,j))*(N(j,7))^2); %Calcula valores de la diagonal de la matriz L
                end
    lt=lt+1;

                if j~=t
                    JA(F,Q)=N(j,7)*N(t,7)*(real(Ybarra(j,t))*sin(N(j,8)-N(t,8))-imag(Ybarra(j,t))*cos(angle(N(j,8)-N(t,8)))); %Calcula valores de la triangula superior y inferior de la matriz L 

                end
                    Q=1+Q;
                    
             end
            end
           F=F+1;
             
        end 
       
end
%_______________________________

for j=1:FN
    if N(j,2)==2
         for t=1:FN
                                                                                             %Faltan restricciones
         P(i,1)=(-1)*(N(j,7)*N(t,7)*abs(Ybarra(j,t))*sin(N(t,8)-N(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
         Pn(i,1)=P(i,1)+Pn(i,1);
            
         end   
         i=i+1;
    end

end

deltaP=Po-Pn;
VA=VAN;
V0=V0N;


i=1;
deltasV=(inv(JA)*deltaP);
for j=1:FN
    if N(j,2)==2||N(j,2)==1
        Vant(i,1)=N(j,8);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        N(j,8)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
end

for j=1:FN
    if N(j,2)==2
        Vant(i,1)=N(j,7);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        N(j,7)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
    
end
if max(abs(deltasV))<error
    fprintf ("El problema converge en la iteraciòn # %i\n",X);
      
    for j=1:FN
    if N(j,2)==0
       for b=1:FN  
  Ns=(-1)*(N(j,7)*N(b,7)*abs(Ybarra(j,b))*sin(N(b,8)-N(j,8)+angle(Ybarra(j,b))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qs=Ns+Qs;
       end
    end

    end
    break
end

end
fprintf("Solucion por el mètodo de NR clasico\n");
% Ybarra %% Ybarra en polares con la estructura (magnitud)+(angulo)j
% Po
% Pn
% deltaP
Qcalculadas
deltasV
if max(abs(deltasV))>error

fprintf("el problema no converge, con ese nùmero de iteraciones, las respuesta mostrada se encuentra en la iteración %i\n", X);
end
Qs 
Vfinal=deltasV+Vant
% Vant
% N
V0=wy(:,7);
VA=wy(:,8);
VAN=VA;
V0N=V0;
Ns=0;
Qs=0;
for X=1:iterMAX
%Obteniendo la potencia activa y reactiva
FP=0;        % # Filas de la matriz de potencias calculadas
NPQ=0;       % # De nodos PQ
NPV=0;       % # De nodos PV
for j=1 :FN  % Se crea la matriz de potencias y su dimencion depende de el # de nodos pv y pq
    H=wy(j,2);FP=H+FP;
    if H==2
        NPQ=NPQ+1;
    end
    if H==1
        NPV=NPV+1;
    end
end
magnitud=(NPQ+NPV);
P=zeros(FP,1);
Pn=zeros(FP,1);
Po=zeros(FP,1);
JA=zeros(FP,FP);
Vant=zeros(FP,1);%% usada para capturar los valores angulo y tencion en los nodos PV y PQ de la matriz wy

Pcalculadas=zeros(magnitud,1);
Qcalculadas=zeros(magnitud,1);
%Creando matrizes para almacenar los valores actualizados de magnitud y
%angulo de tension
%CREANDO VECTOR DE POTENCIAS INICIALES
V=1; 
for j=1:FN
    if wy(j,2)==2||wy(j,2)==1
        Po(V,1)=wy(j,3)-wy(j,5);
        V=V+1;
        
    end
end
for j=1:FN
    if wy(j,2)==2
        Po(V,1)=wy(j,4)-wy(j,6);
        V=V+1;
    end
end
%_______________________________
i=1;% Contador para llevar la cuenta de los numeros que colocan en la matriz P
F=1; %Controla la posicion de las filas del jacobiano
m=1;

for j=1:FN

        if wy(j,2)==1 || wy(j,2)==2 
       for t=1:FN  
  P(m,1)=(-1)*(wy(j,7)*wy(t,7)*abs(Ybarra(j,t))*sin(wy(t,8)-wy(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qcalculadas(m,1)=P(m,1)+Qcalculadas(m,1);
       end
       
       m=m+1;
        end

      if wy(j,2)==1 || wy(j,2)==2 
  for t=1:FN
       P(i,1)=wy(j,7)*wy(t,7)*abs(Ybarra(j,t))*cos(wy(t,8)-wy(j,8)+angle(Ybarra(j,t)));
  Pcalculadas(i,1)=P(i,1)+Pcalculadas(i,1);
  end
      end
      
    if wy(j,2)==2||wy(j,2)==1
        for t=1:FN
        P(i,1)=wy(j,7)*wy(t,7)*abs(Ybarra(j,t))*cos(wy(t,8)-wy(j,8)+angle(Ybarra(j,t)));% hallando P calculada. El valor de p calculada queda en Pn y no en P
        Pn(i,1)=P(i,1)+Pn(i,1);
        end
        i=i+1;
    end
Q=1;   %Controla la posicion de las columnas del jacobiano
if wy(j,2)==2||wy(j,2)==1
    for t=1:FN   
        if wy(t,2)==2||wy(t,2)==1
            if j==t
                JA(F,Q)=-imag(Ybarra(j,j)); %Calcula valores de la diagonal de la matriz H 
            end
            if j~=t
               JA(F,Q)=-imag(Ybarra(j,t)); %Calcula valores de la triangula superior y inferior de la matriz H 
            end                
        Q=1+Q;
         end
    end
    for t=1:FN 
        if wy(t,2)==2
            if j==t
                JA(F,Q)=0; %Calcula valores de la diagonal de la matriz N crow
            end
            if j~=t
                JA(F,Q)=0; %Calcula valores de la triangula superior y inferior de la matriz N crow
            end 
            Q=1+Q;
        end
    end
    F=F+1;
end 

end
%_______________________________

for j=1:FN
    Q=1;   %Controla la posicion de las columnas del jacobiano

    if wy(j,2)==2
        o=1;
    for t=1:FN   
        if wy(t,2)==2||wy(t,2)==1

            if j==t
                JA(F,Q)=0; %Calcula valores de la diagonal de la matriz J crow
            end
              o=o+1;
            if j~=t
                JA(F,Q)=0; %Calcula valores de la triangula superior y inferior de la matriz J crow
            end
        Q=1+Q;

        end

    end
    lt=1;
        for t=1:FN 
            if wy(t,2)==2
            if j==t
                JA(F,Q)=-imag(Ybarra(j,j)); %Calcula valores de la diagonal de la matriz L
            end
            lt=lt+1;
            if j~=t
                JA(F,Q)=-imag(Ybarra(j,t)); %Calcula valores de la triangula superior y inferior de la matriz L 

            end

Q=1+Q;
                    
end
end
F=F+1;
             
end 

end
%_______________________________

for j=1:FN
    if wy(j,2)==2
         for t=1:FN
                                                                                             %Faltan restricciones
         P(i,1)=(-1)*(wy(j,7)*wy(t,7)*abs(Ybarra(j,t))*sin(wy(t,8)-wy(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
         Pn(i,1)=P(i,1)+Pn(i,1);
            
         end   
         i=i+1;
    end

end

deltaP=Po-Pn;
VA=VAN;
V0=V0N;

i=1;
deltasV=(inv(JA)*deltaP);
for j=1:FN
    if wy(j,2)==2||wy(j,2)==1
        Vant(i,1)=wy(j,8);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        wy(j,8)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
end

for j=1:FN
    if wy(j,2)==2
        Vant(i,1)=wy(j,7);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        wy(j,7)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
    
end
if max(abs(deltasV))<error
    fprintf("Solucion por el mètodo de NR desacoplado rapido\n");
    fprintf ("El problema converge en la iteraciòn # %i\n ",X);
      
    for j=1:FN
    if wy(j,2)==0
       for b=1:FN  
  Ns=(-1)*(wy(j,7)*wy(b,7)*abs(Ybarra(j,b))*sin(wy(b,8)-wy(j,8)+angle(Ybarra(j,b))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qs=Ns+Qs;
       end
    end
    end
    break
end

end


% Ybarra %% Ybarra en polares con la estructura (magnitud)+(angulo)j
% Po
% Pn
% deltaP

Qcalculadas
deltasV
if max(abs(deltasV))>error
fprintf("Solucion por el mètodo de NR desacoplado rapido\n")
fprintf("el problema no converge, con ese nùmero de iteraciones, las respuesta mostrada se encuentra en la iteración %i\n", X);
end
Qs 

% Vant
% wy
%crow

VFinish=Vant+deltasV
V0=U(:,7);
VA=U(:,8);
VAN=VA;
V0N=V0;
Ns=0;
Qs=0;
for X1=1:iterMAX
%Obteniendo la potencia activa y reactiva
FP=0;        % # Filas de la matriz de potencias calculadas
NPQ=0;       % # De nodos PQ
NPV=0;       % # De nodos PV
for j=1 :FN  % Se crea la matriz de potencias y su dimencion depende de el # de nodos pv y pq
    H=U(j,2);FP=H+FP;
    if H==2
        NPQ=NPQ+1;
    end
    if H==1
        NPV=NPV+1;
    end
end
magnitud=(NPQ+NPV);
P=zeros(FP,1);
Pn=zeros(FP,1);
Po=zeros(FP,1);
JA=zeros(FP,FP);
Vant=zeros(FP,1);%% usada para capturar los valores angulo y tencion en los nodos PV y PQ de la matriz U

Pcalculadas=zeros(magnitud,1);
Qcalculadas=zeros(magnitud,1);
%Creando matrizes para almacenar los valores actualizados de magnitud y
%angulo de tension
%CREANDO VECTOR DE POTENCIAS INICIALES
V=1; 
for j=1:FN
    if U(j,2)==2||U(j,2)==1
        Po(V,1)=U(j,3)-U(j,5);
        V=V+1;
        
    end
end
for j=1:FN
    if U(j,2)==2
        Po(V,1)=U(j,4)-U(j,6);
        V=V+1;
    end
end
%___________
i=1;% Contador para llevar la cuenta de los numeros que colocan en la matriz P
F=1; %Controla la posicion de las filas del jacobiano
m=1;

for j=1:FN

        if U(j,2)==1 || U(j,2)==2 
       for t=1:FN  
  P(m,1)=(-1)*(U(j,7)*U(t,7)*abs(Ybarra(j,t))*sin(U(t,8)-U(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qcalculadas(m,1)=P(m,1)+Qcalculadas(m,1);
       end
       
       m=m+1;
        end

      if U(j,2)==1 || U(j,2)==2 
  for t=1:FN
       P(i,1)=U(j,7)*U(t,7)*abs(Ybarra(j,t))*cos(U(t,8)-U(j,8)+angle(Ybarra(j,t)));
  Pcalculadas(i,1)=P(i,1)+Pcalculadas(i,1);
  end
      
      end
      
    if U(j,2)==2||U(j,2)==1
        for t=1:FN
        P(i,1)=U(j,7)*U(t,7)*abs(Ybarra(j,t))*cos(U(t,8)-U(j,8)+angle(Ybarra(j,t)));% hallando P calculada. El valor de p calculada queda en Pn y no en P
        Pn(i,1)=P(i,1)+Pn(i,1);
        end
        i=i+1;
    end
    
        Q=1;   %Controla la posicion de las columnas del jacobiano
        
        if U(j,2)==2||U(j,2)==1
        
        for t=1:FN   
            if U(t,2)==2||U(t,2)==1
            
                if j==t
                
                    JA(F,Q)=-Qcalculadas(F,1)-(imag(Ybarra(j,j))*(U(j,7))^2); %Calcula valores de la diagonal de la matriz H 
                    
                end
                
                if j~=t
                
                   JA(F,Q)=U(j,7)*U(t,7)*(real(Ybarra(j,t))*sin(U(j,8)-U(t,8))-imag(Ybarra(j,t))*cos(U(j,8)-U(t,8))); %Calcula valores de la triangula superior y inferior de la matriz H 
                   
                end                
            Q=1+Q;
             end
            end
            
            for t=1:FN 
            
                if U(t,2)==2
                
                if j==t
                
                    JA(F,Q)=0; %Calcula valores de la diagonal de la matriz N crow
                    
                end
                
                if j~=t
                
                    JA(F,Q)=0; %Calcula valores de la triangula superior y inferior de la matriz N crow

                end 
                    Q=1+Q;
             end
            end
           F=F+1;
        end 
        
end

%___________

for j=1:FN
        Q=1;   %Controla la posicion de las columnas del jacobiano
        
        if N(j,2)==2
            o=1;
        for t=1:FN   
            if N(t,2)==2||N(t,2)==1

                
                if j==t
                    JA(F,Q)=0; %Calcula valores de la diagonal de la matriz J crow
                end
                  o=o+1;
                if j~=t
                    JA(F,Q)=0; %Calcula valores de la triangula superior y inferior de la matriz J crow
                end
            Q=1+Q;
         
            end
            
        end

        lt=1;

            for t=1:FN 
                if N(t,2)==2
                if j==t
                    JA(F,Q)=Qcalculadas(lt,1)-(imag(Ybarra(j,j))*(U(j,7))^2); %Calcula valores de la diagonal de la matriz L
                end
    lt=lt+1;

                if j~=t
                    JA(F,Q)=U(j,7)*U(t,7)*(real(Ybarra(j,t))*sin(U(j,8)-U(t,8))-imag(Ybarra(j,t))*cos(angle(U(j,8)-U(t,8)))); %Calcula valores de la triangula superior y inferior de la matriz L 

                end
                    Q=1+Q;
                    
             end
            end
           F=F+1;
             
        end 
       
end
%___________


for j=1:FN
    if U(j,2)==2
         for t=1:FN
                                                                                             %Faltan restricciones
         P(i,1)=(-1)*(U(j,7)*U(t,7)*abs(Ybarra(j,t))*sin(U(t,8)-U(j,8)+angle(Ybarra(j,t))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
         Pn(i,1)=P(i,1)+Pn(i,1);
            
         end   
         i=i+1;
    end

end

deltaP=Po-Pn;
VA=VAN;
V0=V0N;


i=1;
deltasV=(inv(JA)*deltaP);
for j=1:FN
    if U(j,2)==2||U(j,2)==1
        Vant(i,1)=U(j,8);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        U(j,8)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
end

for j=1:FN
    if U(j,2)==2
        Vant(i,1)=U(j,7);% hallando P calculada. El valor de p calculada queda en Pn y no en P
        U(j,7)=Vant(i,1)+deltasV(i,1);
        i=i+1;
    end
    
end
if max(abs(deltasV))<error
    fprintf("Método de NR desacoplado\n");
    fprintf ("El problema converge en la iteraciòn # %i\n",X1);
      
    for j=1:FN
    if U(j,2)==0
       for b=1:FN  
  Ns=(-1)*(U(j,7)*U(b,7)*abs(Ybarra(j,b))*sin(U(b,8)-U(j,8)+angle(Ybarra(j,b))));%Hallando Q calculada. El valor de Q calculada queda en Pn y no en P
  Qs=Ns+Qs;
       end
    end

    end
    break
end

end
% Ybarra %% Ybarra en polares con la estructura (magnitud)+(angulo)j
% Po
% Pn
% deltaP

Qcalculadas
deltasV
if max(abs(deltasV))>error
   fprintf("Método de NR desacoplado\n");
fprintf("el problema no converge, con ese nùmero de iteraciones, las respuesta mostrada se encuentra en la iteración %i\n", X1);
end
Qs 
VFinish=Vant+deltasV
% Vant
% U