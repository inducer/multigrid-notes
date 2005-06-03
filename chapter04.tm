<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|3><section|The W-Cycle>

  Let <with|mode|math|V<rsub|0>\<subset\>V<rsub|1>\<subset\>\<cdots\>\<subset\>V<rsub|J>>
  be a nested sequence of finite element spaces.

  <center|<tabular|<tformat|<table|<row|<cell|<with|mode|math|Q<rsub|j>:V<rsub|J>\<rightarrow\>V<rsub|j>>>|<cell|<with|mode|math|L<rsup|2>>
  projection>|<cell|>>|<row|<cell|<with|mode|math|P<rsub|j>:V<rsub|J>\<rightarrow\>V<rsub|j>>>|<cell|Galerkin
  projection>|<cell|<with|mode|math|P<rsub|j>=A<rsup|-1><rsub|j>Q<rsub|j>A<rsub|j+1>>>>|<row|<cell|<with|mode|math|A<rsub|j>:V<rsub|j>\<rightarrow\>V<rsub|j>>>|<cell|Galerkin
  approximation of the operator>|<cell|<with|mode|math|<ip|A<rsub|j>v<rsub|j>|w<rsub|j>><rsub|0>=a(v<rsub|j>,w<rsub|j>)>>>|<row|<cell|<with|mode|math|R<rsub|j>:V<rsub|j>\<rightarrow\>V<rsub|j>>>|<cell|The
  smoother>|<cell|<with|mode|math|Q<rsub|j>A<rsub|j+1>v<rsub|j>=A<rsub|j>v<rsub|j><with|color|red|?????>>>>>>>>

  We define the W-cycle preconditioner <with|mode|math|B<rsub|j>:V<rsub|j>\<rightarrow\>V<rsub|j>>
  recursively by <with|mode|math|B<rsub|0>=A<rsup|-1><rsub|0>> and
  <with|mode|math|c<rsub|j>=B<rsub|j>r<rsub|j>> (using the residual
  <with|mode|math|r<rsub|j>>) is given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|c<rsub|j><rsup|0>>|<cell|=>|<cell|0>>|<row|<cell|c<rsub|j><rsup|i>>|<cell|=>|<cell|c<rsub|j><rsup|i-1>+R<rsub|j>(r<rsub|j>-A<rsub|j>c<rsub|j><rsup|i-1>)<space|1em>i=1,\<ldots\>,m>>|<row|<cell|c<rsub|j><rsup|m+1>>|<cell|=>|<cell|c<rsub|j><rsup|m>B<rsub|j-1>Q<rsub|j-1>(r<rsub|j>-A<rsub|j>c<rsub|j><rsup|m>)>>|<row|<cell|c<rsub|j><rsup|m+2>>|<cell|=>|<cell|c<rsub|j><rsup|m+1>B<rsub|j-1>Q<rsub|j-1>(r<rsub|j>-A<rsub|j>c<rsub|j><rsup|m+1>)>>|<row|<cell|c<rsub|j>>|<cell|\<assign\>>|<cell|c<rsub|j><rsup|m+2>.>>>>
  </eqnarray*>

  Consider the linear iteration

  <\equation*>
    u<rsub|J><rsup|k>=u<rsub|J><rsup|k-1>+B<rsub|J>(f<rsub|J>-A<rsub|J>u<rsub|J><rsup|k-1>)
  </equation*>

  for a start vector <with|mode|math|u<rsub|J><rsup|0>\<in\>V<rsub|J>>.

  <big-figure|W-Zykel<with|color|red|Bild1>|>

  <subsection|Implementation>

  Start with <with|mode|math|\<b-u\><rsub|J>\<in\>\<bbb-R\><rsup|N<rsub|J>>>,
  compute <with|mode|math|\<b-r\><rsup|J>=\<b-f\><rsub|J>-\<b-A\><rsub|J>\<b-u\><rsub|J>>.
  Then, linear solver for <with|mode|math|k=1,2,3,\<ldots\>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-c\><rsub|J>>|<cell|=>|<cell|\<b-B\><rsub|J>\<b-r\><rsub|J>,>>|<row|<cell|\<b-u\><rsub|J>>|<cell|\<assign\>>|<cell|\<b-u\><rsub|J>+\<b-c\><rsub|J>,>>|<row|<cell|\<b-r\><rsub|J>>|<cell|\<assign\>>|<cell|\<b-r\><rsub|J>-\<b-A\><rsub|J>\<b-c\><rsub|J>,<space|1em><with|mode|text|(computed
    together with <with|mode|math|\<b-c\><rsub|J>>)>>>>>
  </eqnarray*>

  where <with|mode|math|\<b-B\><rsub|0>=\<b-A\><rsup|-1><rsub|0>> and
  <with|mode|math|\<b-B\><rsub|J>> ist defined recursively by, for
  <with|mode|math|j\<in\>{1,\<ldots\>,J}> and given
  <with|mode|math|\<b-r\><rsub|J>>, start with
  <with|mode|math|\<b-c\><rsub|J>=0>:\ 

  for <with|mode|math|i=1,\<ldots\>,m>:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-v\><rsub|j>>|<cell|\<assign\>>|<cell|\<b-R\><rsub|j>\<b-r\><rsub|j>>>|<row|<cell|\<b-c\><rsub|j>>|<cell|\<assign\>>|<cell|\<b-c\><rsub|j>+\<b-v\><rsub|j>>>|<row|<cell|\<b-r\><rsub|j>>|<cell|\<assign\>>|<cell|\<b-r\><rsub|j>-\<b-A\><rsub|j>\<b-v\><rsub|j>>>>>
  </eqnarray*>

  for <with|mode|math|i=1,\<ldots\>,\<gamma\>>:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-r\><rsub|j-1>>|<cell|=>|<cell|\<b-I\><rsub|j><rsup|T>\<b-r\><rsub|j>>>|<row|<cell|\<b-v\><rsub|j-1>>|<cell|=>|<cell|\<b-B\><rsub|j-1>\<b-r\><rsub|j-1>>>|<row|<cell|\<b-v\><rsub|j>>|<cell|=>|<cell|\<b-I\><rsub|j>\<b-v\><rsub|j-1>>>|<row|<cell|\<b-c\><rsub|j>>|<cell|\<assign\>>|<cell|\<b-c\><rsub|j>+\<b-v\><rsub|j>>>|<row|<cell|\<b-r\><rsub|j>>|<cell|\<assign\>>|<cell|\<b-r\><rsub|j>-\<b-A\><rsub|j>\<b-v\><rsub|j>>>>>
  </eqnarray*>

  Error propagation: (observe that <with|mode|math|\<gamma\>=2> for the
  ``typical'' W-cycle)

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-e\><rsup|k+1><rsub|J>>|<cell|=>|<cell|(id-\<b-B\><rsub|J>\<b-A\><rsub|J>)\<b-e\><rsup|k><rsub|J>>>|<row|<cell|>|<cell|=>|<cell|(id-\<b-I\><rsub|j>\<b-B\><rsub|j-1>\<b-I\><rsub|j><rsup|T>\<b-A\><rsub|j>)<rsup|2>(id-\<b-R\><rsub|j>\<b-A\><rsub|j>)<rsup|m>\<b-e\><rsup|k><rsub|J>>>>>
  </eqnarray*>

  For convenience reasons, we now look for
  <with|mode|math|e<rsub|J>\<in\>V<rsub|J>> instead of using the vector
  representation

  <\eqnarray*>
    <tformat|<table|<row|<cell|e<rsup|k+1>>|<cell|=>|<cell|(id-B<rsub|J>A<rsub|J>)e<rsup|k>>>|<row|<cell|>|<cell|=>|<cell|(id-B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)<rsup|2>(id-R<rsub|j>A<rsub|j>)<rsup|m>>>|<row|<cell|<with|mode|text|special
    case: two-level, <with|mode|math|j=1>:
    >>|<cell|=>|<cell|(id-<wide*|A<rsup|-1><rsub|0>Q<rsub|0>A<rsub|1>|\<wide-underbrace\>><rsub|P<rsub|0>>)>>>>
  </eqnarray*>

  <with|mode|math|N<rsub|j>> is the number of unknows at the
  <with|mode|math|j>th level. The operation count for computing
  <with|mode|math|\<b-c\><rsub|j>=\<b-B\><rsub|J>\<b-r\><rsub|J>> is

  <\equation*>
    O(B<rsub|j>)=m*C*N<rsub|j>+\<gamma\>(C*N<rsub|j>+O(B<rsub|j-1>)).
  </equation*>

  Altogether, we have for <with|mode|math|\<gamma\>=2>

  <\equation*>
    <big|sum><rsub|j=0><rsup|J>m*C*<wide*|N<rsub|j>|\<wide-underbrace\>><rsub|O(4<rsup|j>)>+<wide*|\<gamma\><rsup|J+1-j>|\<wide-underbrace\>><rsub|2<rsup|J+1-j>><wide*|O(B<rsub|j?>)|\<wide-underbrace\>><rsub|4<rsup|j>>=C<big|sum><rsub|j=0><rsup|J>N<rsub|J><left|(><frac|1|2><right|)><rsup|j>=2C*N<rsub|J>.
  </equation*>

  <em|Result:> A W-cycle requires <with|mode|math|O(m*N<rsub|J>)> operations.

  <\lemma>
    Let <with|mode|math|u<rsub|J><rsup|0>> be the start iterate and let
    <with|mode|math|u<rsub|J>=A<rsub|J><rsup|-1>f<rsub|J>> be the exact
    solution on the <with|mode|math|J>th level. Then, we have for the error

    <\equation*>
      e<rsub|J><rsup|k>=u<rsub|J><rsup|k>-u<rsub|J>
    </equation*>

    of the W-cycle the <with|mode|math|m> presmoothing steps

    <\equation*>
      e<rsup|k+1><rsub|J>=(id-B<rsub|J>A<rsub|J>)e<rsup|k><rsub|J>
    </equation*>

    with

    <\equation*>
      (id-B<rsub|j>A<rsub|j>)=(id-(2id-B<rsub|j-1>A<rsub|j-1>)B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)(id-R<rsub|j>A<rsub|j>)<rsup|m>
    </equation*>

    for <with|mode|math|j=2,\<ldots\>,J> and
    <with|mode|math|(id-B<rsub|1>A<rsub|1>)=(id-P<rsub|0>)>.
  </lemma>

  <\proof>
    \;

    <\eqnarray*>
      <tformat|<table|<row|<cell|(id-B<rsub|j>A<rsub|j>)>|<cell|=>|<cell|(id-B<rsub|j-1>Q<rsub|j1>A<rsub|j>)<rsup|2>(id-R<rsub|j>A<rsub|j>)>>>>
    </eqnarray*>

    where

    <\eqnarray*>
      <tformat|<table|<row|<cell|(id-B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)<rsup|2>>|<cell|=>|<cell|id-2B<rsub|j-1>Q<rsub|j-1>A<rsub|j>+B<rsub|j-1>Q<rsub|j-1>A<rsub|j>B<rsub|j-1>Q<rsub|j-1>A<rsub|j>>>|<row|<cell|>|<cell|=>|<cell|id-<wide*|(2id-B<rsub|j-1>A<rsub|j-1>)|\<wide-underbrace\>><rsub|id+(id+B<rsub|j-1>A<rsub|j-1>)>B<rsub|j-1>Q<rsub|j-1>A<rsub|j>>>>>
    </eqnarray*>

    where we used the auxiliary caclulation

    <\eqnarray*>
      <tformat|<table|<row|<cell|<ip|Q<rsub|j-1>A<rsub|j>v<rsub|j-1>|w<rsub|j>||>>|<cell|=>|<cell|<ip|A<rsub|j>v<rsub|j-1>|Q<rsub|j-1>w<rsub|j>||>=a(v<rsub|j-1>,Q<rsub|j-1>w<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|<ip|A<rsub|j-1>v<rsub|j-1>|Q<rsub|j-1>w<rsub|j>||>=<ip|A<rsub|j-1>v<rsub|j-1>|w<rsub|j>||>.>>>>
    </eqnarray*>
  </proof>

  <\theorem>
    Let <with|mode|math|<l2norm|(id-P<rsub|j-1>)(id-R<rsub|j>A<rsub|j>)<rsup|m>|0|>\<leqslant\>C<rsub|m>>
    be small enough (yet independent of <with|mode|math|j>) and let
    <with|mode|math|<l2norm|(id-R<rsub|j>A<rsub|j>)<rsup|m>|0|>\<leqslant\>C<rsub|S>>
    (independent of <with|mode|math|m>). Then, we have
    <with|mode|math|<l2norm|id-B<rsub|J>A<rsub|J>||>\<leqslant\>C\<less\>1>
    (independent of <with|mode|math|J>).
  </theorem>

  <\proof>
    <with|mode|math|\<rho\><rsub|j>\<assign\><l2norm|id-B<rsub|j>A<rsub|j>|0|>>,
    we have <with|mode|math|\<rho\><rsub|1>\<leqslant\>C<rsub|m>> small
    enough. Define

    <\equation*>
      K<rsub|j>\<assign\>(id-R<rsub|j>A<rsub|j>).
    </equation*>

    Now, consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|P<rsub|j-1>K<rsub|j><rsup|m>|0|>>|<cell|=>|<cell|<l2norm|(id-P<rsub|j-1>)K<rsub|j><rsup|m>-K<rsub|j><rsup|m>|0|>\<leqslant\>C<rsub|m>+
      C<rsub|S>.>>>>
    </eqnarray*>

    Then

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|id-(2id-B<rsub|j-1>A<rsub|j-1>)B<rsub|j-1>Q<rsub|j-1>A<rsub|j>>>|<row|<cell|>|<cell|=>|<cell|id-P<rsub|J-1>+A<rsup|-1><rsub|j-1>Q<rsub|j-1>A<rsub|j>-(2id-B<rsub|j-1>A<rsub|j-1>)B<rsub|j-1>A<rsub|j-1><wide*|A<rsup|-1><rsub|j-1>Q<rsub|j-1>A<rsub|j-1>|\<wide-underbrace\>><rsub|P<rsub|j-1>>>>|<row|<cell|>|<cell|=>|<cell|(id-P<rsub|j>)+<wide*|(id-2B<rsub|j-1>A<rsub|j-1>+B<rsub|j-1>A<rsub|j>B<rsub|j-1>A<rsub|j-1>)|\<wide-underbrace\>><rsub|(id-B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)<rsup|2>>P<rsub|j-1>.>>>>
    </eqnarray*>

    Now,

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<rho\><rsub|j>>|<cell|=>|<cell|<l2norm|(id-P<rsub|j-1>)K<rsub|j>+(id-B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)<rsup|2><rsub|j-1>K<rsub|j><rsup|m>||>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|m>+\<rho\><rsub|j-1><rsup|2>(C<rsub|m>+C<rsub|S>)>>>>
    </eqnarray*>

    For <with|mode|math|C<rsub|S>=1>, choose <with|mode|math|m> such that
    <with|mode|math|C<rsub|m>\<leqslant\>1/5>. Then,
    <with|mode|math|\<rho\><rsub|j>\<leqslant\>5/3C<rsub|m>\<leqslant\>1/3>
    by induction: <with|mode|math|j=1> is obviously true. Then,
    <with|mode|math|j-1\<rightarrow\>j> works like so:

    <\equation*>
      \<rho\><rsub|j>\<leqslant\>C<rsub|m>+\<rho\><rsub|j-1><rsup|2>(C<rsub|m>+1)\<leqslant\>C<rsub|m>+<frac|5|3>C<rsub|m><frac|1|3><left|(><frac|1|5>+1<right|)>\<leqslant\><frac|5|3>C<rsub|m>\<leqslant\><frac|1|3>.
    </equation*>

    \;
  </proof>

  <\remark*>
    <\enumerate-alpha>
      <item>In practice, always use <with|mode|math|\<gamma\>=1>. In
      practcice, use <with|mode|math|B<rsub|J>> as a preconditioner.

      <item>Variable V-cycle: <with|mode|math|m<rsub|j>=2<rsup|J-j+1>>
      smooting:

      <tabular|<tformat|<table|<row|<cell|1>|<cell|smoothings>|<cell|<with|mode|math|J>>>|<row|<cell|2>|<cell|>|<cell|<with|mode|math|J-1>>>|<row|<cell|4>|<cell|>|<cell|<with|mode|math|J-2>>>|<row|<cell|<with|mode|math|\<vdots\>>>|<cell|>|<cell|<with|mode|math|\<vdots\>>>>>>>

      So, you have the same number of smoothings as W-cycle, and you get a
      bounded condition number <with|mode|math|\<kappa\><rsub|2>(B<rsub|J>A<rsub|J>)>.

      <item>W-cycle result requires only the smoothing property and the
      approximation property. If an <with|mode|math|L<rsup|2>> estimate is
      available (i.e. <with|mode|math|<l2norm|u-u<rsub|h>|0|>\<leqslant\>C*h<rsup|2><l2norm|f|0|>>),
      the approximation property holds automatically. Applications to
      non-conforming or non-nested problems.
    </enumerate-alpha>
  </remark*>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|4|?>>
    <associate|auto-2|<tuple|4.1|?>>
    <associate|auto-3|<tuple|4.1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal||<pageref|auto-2>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>The
      W-Cycle> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>