# Note on vector "intervals" in Hartfiel's *Markov Set-Chains* ch. 2
(c) Marshall Abrams, 8/2017

#### Propositions:

Where [p, q] is an interval in Hartfiel's sense:

     p = x for an x in the interval iff p = q

     q = x for an x in the interval iff p = q

     (In the above case, there is only one x in the interval.)

     If p != q, then (sum p) < 1, and (sum q) > 1.

In other words,

     (sum p) < 1 < (sum q) iff not p = q

#### Proof:

For a vector "interval" (see second page of chapter for definition), the
boundary vectors p and q are stochastic iff p = q, since pi <= qi for
all i.  The only way for p and q to both be stochastic, i.e. sum to 1,
without being the same vector, is if they "cross over".  e.g. (0.25,
0.75) and (0.75, 0.25) are both stochastic, but it's not the case that p
<= q componentwise.  For example, if p is stochastic, and q not equal p,
then q must have at least one element greater or lower than the
corresponding element of p.  It can't be lower, by assumption.  But if
it's greater, and if all of the elements of q are >= elements of p, then
q must sum to more than 1.  So it's not stochastic.

Now, by assumption a Hartfiel interval is never empty.  i.e. there must
be at least one stochastic vector x such that p <= x <= q,
componentwise.  When p = q, x is that same vector.  When p != q, p must
be <= x, componentwise, and q must be >= x, for all (stochastic) x in
the interval.

In fact, if p = x or q = x for some x in the interval, then p = q = x,
and x is the only member of the interval.  For if p is equal some x in
the interval, then any other vectors `x*` in the interval would have all
elements either equal to p's element or > p's element.  However, in that
case `x*` sums to more than 1, so it's not in the interval.

By definition of "interval", each of the x's in an interval of course
sums to 1, and p <= x.  Given what was shown in the preceding paragraph,
if p != q, then for each x in the interval, p has at least one element
that's less than the corresponding element of x.  In that case p sums to
less than 1.  Similarly, for each x, q has at least one element that's
greater than the corresponding element of x, so q sums to more than 1.


