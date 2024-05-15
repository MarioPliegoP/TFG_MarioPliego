%% STATIC CHARACTERIZATION
% Mario Pliego Padilla. 08/02/2024.

%% DESCRIPTION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates 
% GAIN, OFFSET, DNL, INL, MISSING CODES, TRANSFER FUNCTION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% output = sim('untitled.slx');

line = 1.2;
point = 15;

tic
out = sim('ADC_model_jitter_metaestab.slx');
assignin('base','out',out)
toc

%% ERROR WARNINGS

if isempty(out.do.signals.values)==1
        app.errorlabel.Text='Noise too large! Input-Output characteristic cannot be plotted!​';
        return
end

if out.do.signals.values(1)~=0
        app.errorlabel.Text='Noise too large! Input-Output characteristic cannot be plotted!​';
        return
end

for z=1:2^B-1
    brokerep=find(out.do.signals.values==z,1);
    if isempty(brokerep)==1          
       plot(app.transchar,out.ai.signals.values,out.do.signals.values,'b--','LineWidth',line);
       hold(app.transchar,"on")
       plot(app.transchar,AI,DO,'r--','LineWidth',line) 
       xlabel(app.transchar,'AI')
       ylabel(app.transchar,'DO')
       title(app.transchar,'ADC Transfer Characteristic')
       legend(app.transchar,'Simulated Transfer Characteristic','Ideal Transfer Characteristic','Location','southeast')
       hold(app.transchar,"off")
       app.errorlabel.Text='Missing codes detected! Try to lower a1-a3 values!';       
       return 
    end
end

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

for i=1:length(inicio)
    if inicio(i)>final(i)
        app.errorlabel.Text='Noise too large! Input-Output characteristic cannot be plotted!​';
        return
    end
end
for i=1:length(final)-1
     if inicio(i+1)<final(i)
        app.errorlabel.Text= 'Noise too large! Input-Output characteristic cannot be plotted!​';
        return
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


%% REPRESENTATION

plot(app.transchar,out.ai.signals.values,out.do.signals.values,'b--','LineWidth',line);
hold(app.transchar,"on")
plot(app.transchar,AI,DO,'r--','LineWidth',line)
plot(app.transchar,vectorx,vectory,'p-','LineWidth',1)
plot(app.transchar,final,1:2^B-1,'.y','MarkerSize',12)
plot(app.transchar,inicio,0:2^B-2,'.y','MarkerSize',12)
legend(app.transchar,'Simulated Transfer Characteristic','Ideal Transfer Characteristic','Correction: Linking points','Location','southeast')
xlabel(app.transchar,'AI')
ylabel(app.transchar,'DO')
title(app.transchar,'ADC Transfer Characteristic')
hold(app.transchar,"off")

%% STATIC PARAMETERS %%

%% First, we have to calculate T[k]:

T = zeros(1,2^B-1);

for k=1:length(T)
    T(k)=(final(k)+inicio(k))/2;
end

T_ideal = (1:2^B-1)*Q;




%% GAIN AND OFFSET

G = Q*(2^B-2)/(T(2^B-1)-T(1));      % Gain 
app.Gainsol.Text= num2str(abs(G-1)*100,'%.2f');

T1 = T_ideal(1);

Voff = T1 - G*T(1);                          % Offset 
app.Offsetsol.Text=num2str(Voff,'%.2e');

%% TRANFER CHARACTERISTIC CORRECTED:

% G · T[k] + Voff + ε[k] = Q(k-1) + T_1;
% ε[k] residual error -> (final[k] - incio[k])/2?
% ¿Siempre va a salir hacia la derecha, no? (ε[k] > 0)

% epsilon = zeros(1,2^B-1);
% for k=1:2^B-1
%     epsilon(k)=(final(k)-inicio(k))/2;
% end

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

plot(app.correctedgraph,AI,DO,'r--','LineWidth',line)
hold(app.correctedgraph,"on")
plot(app.correctedgraph,vectorx2,vectory2,'b--','LineWidth',line)
legend(app.correctedgraph,'Theoretical Transfer Characteristic','Gain&Offset Corrected Characteristic','Location','southeast')
xlabel(app.correctedgraph,'AI')
ylabel(app.correctedgraph, 'DO')
hold(app.correctedgraph,"off")

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
app.INLSol.Text=sprintf('%5.2f ± %5.2f',INL_max,INL_tolerance);


plot(app.INL,0:2^B-1,[0,INL])
title(app.INL,sprintf('INL(k)'))
legend(app.INL,sprintf('|INL|_{max}= (%3.2f ± %3.2f) LSB',INL_max, INL_tolerance))
xlabel(app.INL,'Code')
ylabel(app.INL,'INL(k)')
grid(app.INL,'on')


%% DNL(k)   (k -> Desde 0 hasta 15).  DNL(0) = 0 (no definido). DNL(15) = 0

W=zeros(1,2^B -1);
for k=1:(2^B - 2)
    W(k)= G*(T(k+1) - T(k));
end

DNL = (W - Q)./Q;
DNL(2^B-1)=0;
DNL_max = max(abs(DNL));

for i=1:length(DNL)
    if abs(DNL(i)) > 0.9 
        fprintf('Missing codes detected! Static ADC performance cannot be determined!​')
        close(figure(1),figure(3),figure(4),figure(5))
        return
    end
end
DNL_tolerance = INL_tolerance;

app.DNLSol.Text=sprintf('%5.2f ± %5.2f',DNL_max,DNL_tolerance);

plot(app.DNL,0:2^B-1,[0,DNL])
title(app.DNL,sprintf('DNL(k)'))
xlabel(app.DNL,'Code')
ylabel(app.DNL,'DNL(k)')
legend(app.DNL,sprintf('|DNL|_{máx}= (%3.2f ± %3.2f) LSB',DNL_max, DNL_tolerance))
grid(app.DNL,'on')
%% -------------------------------------------------- END  -------------------------------------------------- %%