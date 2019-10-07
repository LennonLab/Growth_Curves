A <- 1
umax <- 0.01
L <- 100
b0 <- 0.1
t <- seq(1, 500, 2)

s <- b0+A*exp(-exp(umax*exp(1)*(L-t)/A+1))

plot(t, s, pch = 16, col = "gray", ylim = c(0, 1.2))
abline(v = L, lty = 2, col = "gray", lwd = 2)
abline(h = A, lty = 2, col = "gray", lwd = 2)
abline(h = A + b0, lty = 2, col = "blue", lwd = 2)
abline(-A + b0, umax, lty = 2, col = "red", lwd = 2)
abline(v = L - exp(1)*(b0/umax), lty = 2, col = "blue", lwd = 2)
abline(h = 0, lty = 3, lwd = 2)
