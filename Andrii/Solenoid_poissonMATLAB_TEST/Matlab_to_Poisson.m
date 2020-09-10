stest=fopen('Setup.am','r');

%Einlesen der Datei Setup.am, in der Grundlegeneden Parameter definiert sind und zu scannende Parameter markiert sind. Ich verwende meist die Markierung mit zwei Unterstrichen vor und hinter dem Namen, um die Stellen in der Ursprungsdateien zu zu definieren, die mit Zahlen ersetzten werden sollen.

M=fscanf(stest,'%c'); %Matlabübersetzung der Datei in ein Char-Array

fclose(stest);

%Schreiben der Run-datei auf Grundlage des veränderten char arrays M:

fid=fopen('ScanRun.am','w');

fprintf(fid,'%s',M);

fclose(fid);

% Automesh durchlaufen lassen:

        system('ScanRun.am');

%Poiison ausführen

        system(' poisson SCANRUN.T35');

        pause(5),

%SF7 ausführen und eine Ausgabedatei erzeugen: (in Grid.in7 wird definiert, ob entlang einer line oder über einer Fläche das berechnete Feld ausgeben wird; s. Anhänge)

        system('SF7 ScanRun.T35 Grid.in7');

%% Einlesen der SF7 Datei:  

        fid=fopen('OUTSF7.txt','r');

T=fscanf(fid,'%c');

fclose(fid);

%um den header der OUTSF7.txt Datei zu umgehen hat bei mir folgendes super funktioniert:

[~,a]=regexpi(T,'Index');

out=sscanf(T(a+2:end),'%f');
data = [out(2:9:end),out(4:9:end)]'; % B-Feld, 1te Zeile Z-Achse (cm), 2te B-Feld Z komp. (G)
% Create txt file
fileID = fopen('data.txt','w');
fprintf(fileID,'%e %e\n',data);
fclose(fileID);