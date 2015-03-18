# Introduction #

_From an article written on the old Starscream Wordpress blog at www.jayjaybillings.org/starscream dated 20100925_

I’ve started compiling a list of the things that I like, dislike and learned from my experience developing Starscream 1.x. As I mentioned in my last post, it is has been a couple of years (and now a few months more) since I originally released Starscream and I’ve had the opportunity to learn from some really great people.

# Likes #

**Starscream 1.x is a simple code. I leveraged as many existing packages as possible and I wrote it in C with some minimalist python extensions (never released).**Starscream 1.x was reasonably well documented, given that it originated with a grad school project and evolved into a spare time project.
**I tried to focus on quality during development. So… this was never really made manifest in the code itself because I didn’t know much about software quality in those days, but I tried to provide support, documentation, a test problem, etc.**Starscream 1.x has a very, very, very nice Gadget2 reader. I adapted the routine from one of Volker Springel’s scripts in Gadget2 and it turned out to work pretty well.

# Dislikes #

**Well, numero uno: it’s a crappy code. I didn’t know any better then, but I do now and it eats at me.**Starscream 1.x is not even remotely extensible. It’s a purely procedural code in the strictest sense and is really hard-wired (which I didn’t think was the case at the time…). I started having bad feelings about things right befSore I dropped the initial release when I was playing with the random number generator and thinking about how a user could change out galaxies.
**Starscream 1.x has no tests. There is a “make test” case, but it only checks compilation and linkage, not functionality or validity.**Starscream 1.x only works on Linux. It seemed like a good idea at the time, but then I realized that most of the potential users are actually running Windows. Ooops.

# Lessons Learned #

**I learned the value of utilizing existing tools. That’s something I should continue, to within reason, on version 2.0.**I learned the importance of a good build system… to the extent that Autotools is a good build system (CMake is soooo much better). Building could be a pain unless I had previously tested Starscream on that type of system. Autotools helped a lot with that, although Autotools can be a pain too.
**The python bindings that I never released were awesome. Having the ability to dynamically load and test things was very enabling.**

That’s all for now… dinner time. More as I think of it.