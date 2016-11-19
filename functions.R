
  #### Sensor coefficents

rinko_coefs = list(serial = "#0263 ARO-CAV", A = -42.34162, B = +127.6475, C = -0.3677435, D = +0.01137, E = +0.0046, F = +7.57e-05)
optode_coefs = read.csv("optode_coefs.csv")

  #### Sensor functions

optode.analogtemp <- function(v){
  return( ((v * 45) / 5) - 5 )
}
optode.analogDphase <- function(v){
  return( 10 + (v / 5) * 60 )
}
optode.phaseCalc <- function(DPhase, Temp, coefs){
  # for mkl optodes 3830 & 3835
    with(coefs, {
      print(paste("using foil batch coefs", batch[1]))
      (C0[1]+C0[2]*Temp+C0[3]*Temp^2+C0[4]*Temp^3) +
      (C1[1]+C1[2]*Temp+C1[3]*Temp^2+C1[4]*Temp^3) *
      DPhase+(C2[1]+C2[2]*Temp+C2[3]*Temp^2+C2[4]*Temp^3) *
      DPhase^2+(C3[1]+C3[2]*Temp+C3[3]*Temp^2+C3[4]*Temp^3) *
      DPhase^3+(C4[1]+C4[2]*Temp+C4[3]*Temp^2+C4[4]*Temp^3) *
      DPhase^4
    })
}
par_from_voltage <- function(x, factor, offset){
    return(factor * exp(offset * x))
}
rinko_temp <- function(V, tC = list(A = -5.326887e+00, B = +1.663288e+01, C = -2.123968e+00, D = +4.543014e-01)){
    # RINKO III defaults to #0263 ARO-CAV

  temp = tC$A + tC$B * V + tC$C * V^2 + tC$D * V^3
  return(temp)
}
rinko_o2 <- function(V, t, S, oC = list(A = -4.234162e+01,
                                     B = +1.276475e+02,
                                     C = -3.677435e-01,
                                     D = +1.137000e-02,
                                     E = +4.600000e-03,
                                     F = +7.570000e-05),
                     p = 0.1, G = 0, H = 1){

  # V = output voltage
  # t = tempeture from rinko_temp
  # p = in-situ pressure in decibar
  # G & H = RINKO calibration coefs (alpha and beta)

    # RINKO III #0263 ARO-CAV

  P1 = oC$A / (1 + oC$D * (t - 25) + oC$F * (t - 25)^2)
  P2 = oC$B / (V * (1 + oC$D * (t - 25) + oC$F * (t - 25)^2) + oC$C)
  P = P1 + P2

    # G and H are calibration coefs
  DO = G + H * P
    # pressure correction
  d = p * 0.01 # convert from decibar to MPa
  DO = DO * (1 + oC$E * d)

  # from garcia and gordon
    A0 = 2.00856
    A1 = 3.224
    A2 = 3.99063
    A3 = 4.80299
    A4 = 0.978188
    A5 = 1.71069
    B0 = -0.00624097
    B1 = -0.00693498
    B2 = -0.00690358
    B3 = -0.00429155
    C0 = -3.1168E-07
    Ts = log((298.15-t)/(273.15+t))

    Cstar = exp(A0+(A1*Ts)+(A2*Ts^2)+
                    (A3*Ts^3)+(A4*Ts^4)+(A5*Ts^5)+
                    S*(B0+(B1*Ts)+(B2*Ts^2)+(B3*Ts^3))+
                    (C0*S^2))

    DO = ((Cstar * 44.614 * DO) / 100) * 0.0319988 # mg/l
    DO = (DO/31.9988) * 1000 # mg/l to mmol m-3

  return(DO)
}
