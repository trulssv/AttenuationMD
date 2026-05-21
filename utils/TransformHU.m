function MonoHU = TransformHU(Mono, E)
    MuWa = GetMu('water', E);
    MuAir = GetMu('air', E);
    MonoHU = 1000 * (Mono - 0.1*MuWa) / (0.1*MuWa - 0.1*MuAir);
end