# -*- coding: utf-8 -*-
 2 import numpy as np
 3 import matplotlib.pyplot as plt
 4
 5 num_sims = 5 ### display five runs
 6
 7 t_init = 3
 8 t_end  = 7
 9 N      = 1000 ### Compute 1000 grid points
10 dt     = float(t_end - t_init) / N
11 y_init = 0
12
13 c_theta = 0.7
14 c_mu    = 1.5
15 c_sigma = 0.06
16
17 def mu(y, t):
18     """Implement the Ornstein–Uhlenbeck mu.""" ## = \theta (\mu-Y_t)
19     return c_theta * (c_mu - y)
20
21 def sigma(y, t):
22     """Implement the Ornstein–Uhlenbeck sigma.""" ## = \sigma
23     return c_sigma
24
25 def dW(delta_t):
26     """Sample a random number at each call."""
27     return np.random.normal(loc = 0.0, scale = np.sqrt(delta_t))
28
29 ts    = np.arange(t_init, t_end, dt)
30 ys    = np.zeros(N)
31
32 ys[0] = y_init
33
34 for _ in range(num_sims):
35     for i in range(1, ts.size):
36         t = (i-1) * dt
37         y = ys[i-1]
38         ys[i] = y + mu(y, t) * dt + sigma(y, t) * dW(dt)
39     plt.plot(ts, ys)
40
41 plt.show()
