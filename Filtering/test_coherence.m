% test coherence
clear 

fs = 1;
mtime = 1:1E5;

%signal 1
T1 = [10,150];
a1 = [1,0.2];
p1 = [0 pi/2];
e1 = [0.01 0.02];
f1=0*mtime;
for i=1:length(T1)
   f1=f1+a1(i)*sin(2*pi/T1(i)*mtime+p1(i))+e1(i)*rand(1,length(mtime)); 
end

%signal 2
T2 = [10,150];
a2 = [1,0.2];
p2 = [pi/4 pi/6];
e2 = [0.01 0.02];
f2=0*mtime;
for i=1:length(T1)
   f2=f2+a2(i)*sin(2*pi/T2(i)*mtime+p2(i))+e2(i)*rand(1,length(mtime)); 
end

C = coherence(f1,f2,10000,fs);

p=figure(1),
subplot(2,1,1),semilogx(C.f,C.coherence),ylabel('Coh^2')
subplot(2,1,2),semilogx(C.f,C.phase),ylabel('Phase')
xlabel('f')