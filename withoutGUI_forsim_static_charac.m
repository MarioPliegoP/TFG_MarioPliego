%% Theoretical tranfer characteristic


DO = zeros(2*(2^B),1)';
for k=1:2*(2^B)
    if mod(k,2) == 0
            DO(k)=(k/2)-1;
        else
            DO(k)= (k-1)/2;
    end
end
AI = DO*Q;
AI = AI(2:end);
AI = [AI,(2^B)*Q];


%% Execution of Simulink's model

line = 1.2;
point = 15;

out = sim('ADC_model_jitter_metaestab.slx');

%% CORRECTED SIGN REPRESENTATION

final = zeros(1,2^B-1);
inicio = zeros(1,2^B-1);

for c = 1:length(inicio)
    indice=find(out.do.signals.values==c,1);
    inicio(c)= out.ai.signals.values(indice-1);              
end

for i = 1:length(out.do.signals.values)
    k = out.do.signals.values(i);
    final(k+1) = out.ai.signals.values(i+1);        
    if out.do.signals.values(i+1)==2^B-1
             for h = length(out.do.signals.values):-1:1
               hh = find(out.do.signals.values==2^B-2);
               final(k+1) = out.ai.signals.values(hh(end)+1);         
             end
      break
    end
end


%% 1ST VERSION: JUST LINKING THE POINTS

vectory = zeros((2^B)*2,1);
for k = 1:(2^B)
    vectory(2*k-1) =k-1;
    vectory(2*k)=k-1;
end

vectorx= zeros(1,2*(2^B-1)+2);

for i=2:length(vectorx)-1
    if mod(i,2) == 0
        vectorx(i)=inicio(i/2);
    else
        vectorx(i)=final((i-1)/2);
    end
end

vectorx(1)=out.ai.signals.values(1);
vectorx(length(vectorx))=out.ai.signals.values(end);


%% First, we have to calculate T[k]:

T = zeros(1,2^B-1);

for k=1:length(T)
    T(k)=(final(k)+inicio(k))/2;
end

T_ideal = (1:2^B-1)*Q;


%% GAIN AND OFFSET

G = Q*(2^B-2)/(T(2^B-1)-T(1));      % Gain 
error_G = abs(G-1) * 100;

T1 = T_ideal(1);

Voff = T1 - G*T(1);                          % Offset 

%% TRANFER CHARACTERISTIC CORRECTED:

corrected_AI = G * T + Voff;

vectorx2 = zeros(1,(2^B-1)*2);

for g=1:2:length(vectorx2)
    vectorx2(g)= corrected_AI((g+1)/2);
    vectorx2(g+1) = vectorx2(g);
end

vectory2= zeros(1,(2^B-1)*2+2);

for g=3:length(vectory2)
    if mod(g,2) == 0
        vectory2(g)=(g-2)/2;
    elseif mod(g,2) == 1
        vectory2(g)= (g-1)/2;
    end
end

vectorx2 = [0,vectorx2,FS];

%% INL

INL = ((G*T + Voff - T_ideal)/Q);     % INL
INL_max = max(abs(INL));

INL_tolerance = 0;
for i=1:2^B-1
    if (final(i)-inicio(i))/Q > INL_tolerance
        INL_tolerance = (final(i)-inicio(i))/Q;
    end
end

INL_tolerance = INL_tolerance/2;






