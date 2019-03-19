# HIV Modelling

# Tutorials
  * nice brief tutorial on using LOESS/splines/interpolation as parameters in ODE models https://kingaa.github.io/thid/odes/ODEs_in_R.pdf
    * came from here, maybe more stuff, Aaron King stuff is trustworthy https://kingaa.github.io/thid/
  * something else on spline interpolation https://books.google.com/books?id=ejY8DQAAQBAJ&pg=PA25&lpg=PA25&dq=spline+as+a+parameter+of+ode+model+r+desolve&source=bl&ots=5nJE6LXCvY&sig=BACy3MeiOvy4uwOP5lgOQT8uH9M&hl=en&sa=X&ved=2ahUKEwjchY_6n7LfAhXdwMQHHc6MDY8Q6AEwCHoECAIQAQ#v=onepage&q=spline%20as%20a%20parameter%20of%20ode%20model%20r%20desolve&f=false
  * more on smoothed time-varying parameters: http://www.simecol.de/modecol/

* Model issues:
  * alpha parameter (people can change behavior) - this should probably be added (but let's think about sociologically/culturally what it means...)
  * proportions of pop in each risk category is constant (this could be fixed partially with alpha, but would suffer problems of multicollinearity)
  * not really an issue, but a "gotcha"; *birth* parameter is per-capita, not per-female, also its really the birth at a lag of 15 years multiplied by survival to age 15 (or whatever age of sexual encounters is)
  * all treatment parameters are the same between M/F

* Stuff for future development
  * should probably model partnerships forming/dissolving, need to talk to Sandi and (maybe) Joe
  * demographics should be done more cautiously
  * better mixing matrices (do more research on this)
