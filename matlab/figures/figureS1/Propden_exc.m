function Ptot_exc = Propden_exc(x,PropDens,Pcon,N,Ninh)

Ptot_exc = 0;

for input = (Ninh+1):N
    Ptot_exc = Ptot_exc + PropDens(x, input)*Pcon(input);
end
Ptot_exc = Ptot_exc/sum(Pcon((Ninh+1):N));