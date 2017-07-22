I thought that it might help to make the internal lists at each tick be
lazy.  I thought that maybe this could allow me to run for more
generations out.  This branch was the start of an attemp to do that.  It
doesn't compile yet.

However, I now think that making the internal at-tick lists lazy doesn't
help anything for my current purposes.  When I plot at a tick, I need
all of the contents of this list at once.  There's no reason to make
that lazy.

(Perhaps it could be useful to use a different laziness module that
would allow early generations to disappear once they'd been used.)
