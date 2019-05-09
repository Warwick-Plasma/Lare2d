FUNCTION sound, a

  gamma = 5.0 / 3.0
  cs = SQRT(gamma * a.pressure / a.rho)

  RETURN, cs
END

FUNCTION alfven, a

  va = SQRT((a.bx^2+a.by^2+a.bz^2)/a.rho)

  RETURN, va
END
