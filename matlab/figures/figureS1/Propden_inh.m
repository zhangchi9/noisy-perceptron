function Ptot_inh = Propden_inh(x,PropDens,Pcon,Ninh)

Ptot_inh = 0;

for input = 1:Ninh
    Ptot_inh = Ptot_inh + PropDens(x, input)*Pcon(input);
end
Ptot_inh = Ptot_inh/sum(Pcon(1:Ninh));