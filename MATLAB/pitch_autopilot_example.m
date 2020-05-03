% Pitch Attitude Autopilot
%

G = tf([3], [1 2 5]);
H = tf([1], [1 10]);
GH = G*H;

pause

pzmap(GH)

pause

rlocus(GH)

[Kcr, r] = rlocfind(GH)

Tcr = 2*pi/imag(r(2));

% P control

pause

Kp=0.5*Kcr;

Gc_p=tf([Kp]);

Gcl_p = feedback(series(Gc_p, GH), tf([1],[1]))

% step(Gcl_p)

% PI control

pause

Kp=0.45*Kcr;
Ti=Tcr/1.2;

Gc_pi=tf([Kp*Ti Kp], [Ti 0])

Gcl_pi = feedback(series(Gc_pi, GH), tf([1],[1]))

% step(Gcl_pi)

% PID control

pause

Kp=0.6*Kcr;

Ti=0.5*Tcr;

Td=0.125*Tcr;

Gc_pid= tf([Kp*Td*Ti Kp*Ti Kp], [Ti 0])

Gcl_pid = feedback(series(Gc_pid, GH), tf([1],[1]))

step(Gcl_p, Gcl_pi, Gcl_pid)

GcGH = series(Gc_pid/Kp, GH);

pause

pzmap(GcGH)

pause

rlocus(GcGH)





