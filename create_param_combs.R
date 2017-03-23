# create parameter list
# two_step_binary(ls_size, p, h, w1, w2, f1, f2, fp1_same = True, fp2_same = True)
p <- h <- seq(0.1, 0.9, 0.1)
w1 <- w2 <- c(3, 9, 15)
f1 <- f2 <- c("linear", "exp", "negexp")
r <- 1:100

params <- expand.grid(p, h, w1, w2, f1, f2, r)
names(params) <- c("p", "h", "w1", "w2", "f1", "f2", "r")
write.csv(params, file="two_step_binary/params.csv", row.names = FALSE)

# one_step_binary(ls_size, p, h, w1, f1, fp1_same = True)
params <- expand.grid(p, h, w1, f1, r)
names(params) <- c("p", "h", "w", "f", "r")
write.csv(params, file="one_step_binary/params.csv", row.names = FALSE)
