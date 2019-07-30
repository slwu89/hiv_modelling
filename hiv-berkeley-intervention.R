



int_ranges <- matrix(data = c(
  0.8080,0.8840, # sup1rate_pre2015_int
  0.7000,0.7920, # sup1rate_post2015_int
  0.0028,0.0268, # falloff1_pre2015_int
  0.0438,0.1000 # falloff1_post2015_int
),nrow = 4,ncol = 2,byrow = TRUE,
dimnames = list(c("sup1rate_pre2015_int",
                  "sup1rate_post2015_int",
                  "falloff1_pre2015_int",
                  "falloff1_post2015_int"),
                c("min","max"))
)

