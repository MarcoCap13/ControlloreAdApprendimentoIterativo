
%Tirocinio
%Per prima cosa devo creare un ciclo di controllo dell'impianto
%variabili da implementare: Y_d(reference trajectory) U_k(Input) Y_k(Output)
%e_k(tracking error) = y_d - y_k;


% ---- Traiettoria-->  reference trajectory denotata come y_d sugli appunti
% (pag. 3)
clear all
close all

traiettoria= 0:0.01:10;                                   %è un vettore per identificare la traiettoria che incrementa di 0.01
lunghezza=length(traiettoria);
y_d = zeros(1, lunghezza);
for i=0 :lunghezza-1
    y_d(i+1)=traiettoria(i+1)^3*(4-0.3*traiettoria(i+1))*0.01;
end

figure(1)
plot(traiettoria,y_d);
title('Sistema di riferimento Y_d')

%====================================================
%prova a modificare le matrici A/B che sono quelle che modificano proprio
%lo stato del sistema--> errore input e output modificati
%Variabili che utilizzo per la matrice quadrata A
q=0.8;
s=0.5;
e=1.5;
r=0.01;

%prova1 matrice quadrata
%A=[1 0.1;0.05 1.1];
%prova2 matrice quadrata --> quella che si comporta megio col sistema!
A=[0 1;(q*r/(e*s^2))-1 2-(q*r/(e*s^2))];	 %A è una matrice quadrata
%prova3 matrice quadrata
%A=[0 1;(q*r/e^5)-1 2-(q*r/(e^5))]
B=[0;(r^2)/(e*s^2)];                         %B è un vettore colonna                
%B= [0.45; 0.6];
C=[0 1];
R = 5;                                       %rho --> cambia lo smorzamento
%================================================================
%Studiamo e risolviamo l'equazione di riccati (sist. a tempo discreto) --> eq. differenziale per le
%matrici quadratiche

%legenda:    il nostro K --> è P nel pdf
P(:,:,1001)=zeros(2,2); % P(N)
G(1001,:)=zeros(1,2); % Guadagno
for t=1:1000
    %leggi da pag. 20->22 degli appunti ricordati che K(una matrice) ha
    %soluzioni dell'equazioni di Riccati alle diff. all'indietro.
    %N.B: Q = C'*C
    P(:,:,1001-t)=A'*P(:,:,1001-t+1)*A + C'*C - A'*P(:,:,1001-t+1)*B* ...
        inv(B'*P(:,:,1001-t+1)*B + R)*B'*P(:,:,1001-t+1)*A;
    %il controllo ottimo  in retroazione deve essere lineare il cui
    %guadagno deve essere come l'equazione portata nella pag 22 degli
    %appunti.
    G(1001-t,:) = inv(B'*P(:,:,1001-t)*B+R)*B'*P(:,:,1001-t)*A; % Guadagno di kalman per sistema a  tempo discreto
end
%===============================================================

%Inizializzo le variabili
iterazioni = 10;
lunghezza = 1000;
e_k=zeros(1,1001);                             % Errore
u=zeros(1,1000);                               % Input dell'iterazione precedente
x(:,1001)=zeros(2,1);                          % stato della precedente iterazione
x_k(:,1)=zeros(2,1);                           % stato iniziale corrente

u_k = zeros(iterazioni+1, lunghezza);
u_k(1,:)=ones(1, lunghezza);                   % input corrente
y = zeros(iterazioni, lunghezza);

%===============================================================

%Riprendo  la parte del Pdf riguardo il Tempo Discreto dell'ILC:
%il sistema è formato dalla derivata I Xk(t) = A*Xk(t) + B*Uk(t)
%                                      Yk(t) = C*Xk(t)
%N.B: A,B,C sono sistemi matriciali.
%finchè è richiesta una ripetizione di un certo task, lo stato iniziale
%deve essere resettato per ogni iterazione---> X_k(0) = X_0 per ogni K
%mettiamo il caso di fare 10 iterazioni..

%Creo un ciclo for di controllo per far scendere l'errore a 0 alla fine
%delle 10 iterazioni

for k=1 : iterazioni                                     %-->k = iterazioni
    %Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 1001-by-2.
    %Quest'errore mi è apparso perchè lo spazio vettoriale richiesto è
    %troppo piccolo rispetto al guadagno causato proprio dall'eq. di
    %riccati (1001)
    %creo una nuova matrice di zeri per inserire le iterazioni al suo
    %interno
    indietro(:,1001)=zeros(2,1); % sarebbe l'epsilon(N)
    %prova a scrivere l'equazione alle differenze finite di ordine m 
    for m=1:1000
        %per rappresentare la matrice identica A di ordine n, usa il metodo
        %'eye'
        %risolviamo l'ultima parte dell'eq. di riccati--> RICORDA: di
        %moltiplicare l'errore di tracciamento e_k, se non dovesse
        %funzionare rivediti la parte degli appunti riguardante riccati
        %perchè potrebbero esserci problemi con l'eq. di riccati
        E = B*inv(R)*B';
        indietro(:,1001-m) = inv(eye(2,2) + (P(:,:,1001-m) * E))*(A'* indietro(:,1001-m+1) + C' * e_k(1001-m+1));
    end
    for t=1 : lunghezza                                                                                      
        % ---- Update Law--> falla con u_k e non u_k+1 dove G è la  matrice
        % di guadagno d'apprendimento, u_k è il nostro input, (x_k(:,t)-x(:,t))è
        % l'"innovation term", (1/R)-->[inv(R)] è lo smorzamento del sistema.
        % siccome ho utilizzato l'update law u_k devo prendere l'iterazione
        % precedente e sottrarla al guadagno * l'innovation term
        % -->se non dovesse funzionare riguardati gli appunti a pagina 5 (1.6) e il pdf a pag 24 
        u_k(k,t) = u(t) - (G(t,:) * (x_k(:, t) - x(:,t))) + (inv(R) * B' * indietro(:,t));   %Input  
        x_k(:,t+1) = A*x_k(:,t) + B*u_k(k,t);               % Stato del sistema (vedi pagina 4 degli appunti)
        y(k,t)= C*x_k(:,t);                                 % Output Y_k
        e_k(t)=y_d(t)-y(k,t);                               % Errore di tracciamento --> (traiettoria - output)
        
        %salva gli stati calcolati
        x(:, t) = x_k(:,t);                                 % Salvo lo stato del sistema per la prossima iterazione
        u(t) = u_k(k, t);                                   % Salvo l'input per la prossima iterazione
    end
    x(:,t+1)=x_k(:,t+1);                                    % pongo lo stato della prossima iterazione come iterazione attuale
end



%Nell'output voglio vedere che per ogni iterazione che passa, la triettoria
%n prenderà le 'sembianze' della reference Trajectory
figure(2);
subplot(2,1,1)%Output
plot(traiettoria(1:end-1), y(1,:), traiettoria(1:end-1), y(2,:), traiettoria(1:end-1), y(4,:), traiettoria(1:end-1), y(6,:), traiettoria(1:end-1), y(8,:), traiettoria(1:end-1), y(10,:));
title('Output Y_k per le differenti iterazioni 1-10');
legend('Iterazione 1','Iterazione 2','Iterazione 4','Iterazione 6','Iterazione 8','Iterazione 10');

%Per gli input dovrei aspettarmi che le varie traiettorie tornino a
%0-->ergo  l'inpug U_k è ottimo e non dovremmo migliorarlo nell'iterazione
%successiva.----> FUNZIONA!!
subplot(2,1,2) %Input 
plot(traiettoria(1:end-1), u_k(2,:), traiettoria(1:end-1), u_k(3,:), traiettoria(1:end-1), u_k(4,:), traiettoria(1:end-1), u_k(6,:), traiettoria(1:end-1), u_k(8,:), traiettoria(1:end-1), u_k(10,:));
title('Input U_k per le differenti iterazioni 1-10');
legend('Iterazione 1','Iterazione 2','Iterazione 4','Iterazione 6','Iterazione 8','Iterazione 10');

%Tracciamo l'errore 
figure(3);
plot(traiettoria(1:end-1), e_k(1:end-1));
title('Errore e_k');

%Confornto tra la traiettoria di riferimento con l'ultima iterazione
%Se l'algoritmo dovesse funzionare decentemente, l'ultima iterazione
%dovrebbe riprendere la traiettoria di riferimento.-->FUNZIONA!! le due
%traiettorie sono perfettamente sovrapposte!
figure(4);
plot(traiettoria(1:end-1), y(10,:), traiettoria(1:end-1), y_d(1:end-1));
title('Traiettoria di Riferimento con l ultima iterazione');
legend('Outuput 10','Riferimento');


