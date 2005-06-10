<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|4><section|The symmetric V-cycle>

  Let <with|mode|math|><with|mode|math|V<rsub|0>\<subset\>\<cdots\>\<subset\>V<rsub|J>>
  be a sequence of nested Finite Element Spaces.

  <big-figure|a) <with|color|red|Bild> b) <with|color|red|Bild>|a)... b)...>

  <with|color|red|...Wolfgang...>

  <with|mode|math|K<rsub|j>> is just <em|any> smoother.

  The error propagation is

  <\eqnarray*>
    <tformat|<table|<row|<cell|E<rsub|j>>|<cell|=>|<cell|
    id-B<rsub|j>A<rsub|j>>>|<row|<cell|>|<cell|=>|<cell|K<rsub|j>(id-B<rsub|j-1>Q<rsub|j-1>A<rsub|j>)K<rsub|j><rsup|\<ast\>>>>|<row|<cell|>|<cell|<above|=|!>>|<cell|K<rsub|j>(id-P<rsub|j-1>+E<rsub|j-1>P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>,>>>>
  </eqnarray*>

  where the step marked ``<with|mode|math|!>'' is shown by

  <\eqnarray*>
    <tformat|<table|<row|<cell|B<rsub|j-1>Q<rsub|j-1>A<rsub|j>>|<cell|=>|<cell|P<rsub|j-1>-(A<rsub|j-1><rsup|-1>Q<rsub|j-1>A<rsub|j>)+B<rsub|j-1>A<rsub|j-1>(A<rsub|j-1><rsup|-1>Q<rsub|j-1>A<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|P<rsub|j-1>-<wide*|(id-B<rsub|j>A<rsub|j-1>)|\<wide-underbrace\>><rsub|E<rsub|j-1>>P<rsub|j-1>>>>>
  </eqnarray*>

  For example, at the zeroth level,

  <\eqnarray*>
    <tformat|<table|<row|<cell|E<rsub|0>>|<cell|=>|<cell|0>>|<row|<cell|E<rsub|1>>|<cell|=>|<cell|K<rsub|1>(id-P<rsub|0>)K<rsub|1><rsup|\<ast\>>.>>>>
  </eqnarray*>

  <\theorem>
    Assume the smoothing property

    <\equation>
      <label|eq:vcycle-smoothingprop><enorm|K<rsub|j>v<rsub|j>||2>\<leqslant\>a((id-\<theta\><rsub|j>A<rsub|j>)v<rsub|j>,v<rsub|j>)<space|1em><with|mode|text|with
      <with|mode|math|\<lambda\><rsub|j>\<theta\><rsub|j>=\<omega\>\<in\>(0,1)>>
    </equation>

    and the approximation property

    <\equation>
      <label|eq:vcycle-approxprop><l2norm|v<rsub|j>-P<rsub|j>v<rsub|j>||2>\<leqslant\>C<rsub|p>\<lambda\><rsub|j><rsup|-1><enorm|v<rsub|j>||2>.
    </equation>

    Then, we have for the V-cycle

    <\equation*>
      <enorm|E<rsub|j>||>=<enorm|id-B<rsub|j>A<rsub|j>||>\<leqslant\>\<delta\>\<assign\><frac|C<rsub|p>|C<rsub|p>+\<omega\>>\<less\>1<space|1em>(j=0,1,\<ldots\>J).
    </equation*>
  </theorem>

  <\proof>
    <em|Step 1.> Show by induction

    <\equation*>
      a(E<rsub|j>v<rsub|j>,w<rsub|j>)=a(v<rsub|j>,E<rsub|j>w<rsub|j>).
    </equation*>

    For <with|mode|math|j=0>, <with|mode|math|E<rsub|0>=0><with|mode|math|\<rightarrow\>>bingo.
    Now, for <with|mode|math|j-1\<rightarrow\>j>:

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(E<rsub|j>v<rsub|j>,w<rsub|j>)>|<cell|=>|<cell|a(K<rsub|j>(id-P<rsub|j-1>+E<rsub|j-1>P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>v<rsub|j>,w<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|j-1>+E<rsub|j-1>P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>v<rsub|j>,K<rsub|j><rsup|\<ast\>>w<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|a(K<rsub|j><rsup|\<ast\>>v<rsub|j>,K<rsub|j><rsup|\<ast\>>w<rsub|j>)-a(P<rsub|j-1>K<rsub|j-1><rsup|\<ast\>>v<rsub|j>,<wide*|P<rsub|j-1>|\<wide-underbrace\>><rsub|=id>K<rsub|j-1><rsup|\<ast\>>w<rsub|j>)+<wide*|a(E<rsub|j-1>P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>,<wide*|P<rsub|j-1>|\<wide-underbrace\>><rsub|=id>K<rsup|\<ast\>><rsub|j>w<rsub|j>)|\<wide-underbrace\>><rsub|<with|mode|text|Ind.Ass.:
      >=a(P<rsub|j>,K<rsub|j>v<rsub|j>,E<rsub|j-1>P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>)>>>>>
    </eqnarray*>

    <em|Step 2.> Because of the symmetry of <with|mode|math|E<rsub|j>>,
    <with|mode|math|<enorm|E<rsub|j>||>\<leqslant\>\<delta\>\<Leftrightarrow\>a(E<rsub|j>v<rsub|j>,v<rsub|j>)\<leqslant\>\<delta\><enorm|v<rsub|j>||2>>
    holds. To see this, consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|sup<rsub|v<rsub|j>><frac|a(E<rsub|j>v<rsub|j>,v<rsub|j>)|<enorm|v<rsub|j>||2>>>|<cell|\<leqslant\>>|<cell|sup<rsub|v<rsub|j>,w<rsub|j>><frac|a(E<rsub|j>v<rsub|j>,w<rsub|j>)|<enorm|v<rsub|j>||>*<enorm|v<rsub|j>||>>=sup<rsub|v<rsub|j>><frac|<enorm|E<rsub|j>v<rsub|j>||>|<enorm|v<rsub|j>||>>>>|<row|<cell|>|<cell|=>|<cell|<enorm|E<rsub|j>||>>>|<row|<cell|>|<cell|=>|<cell|sup<rsub|<enorm|v<rsub|j>||>=<enorm|w<rsub|j>||>=1>a(E<rsub|j>v<rsub|j>,w<rsub|j>)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|<enorm|v<rsub|j>||>=<enorm|w<rsub|j>||>=1><frac|1|2>a(E<rsub|j>v<rsub|j>,v<rsub|j>)a(E<rsub|j>w<rsub|j>,w<rsub|j>)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|<enorm|v<rsub|j>||>=1>a(E<rsub|j>v<rsub|j>,v<rsub|j>).>>>>
    </eqnarray*>

    <em|Step 3.> Show <with|mode|math|<l2norm|A<rsub|j>K<rsub|j><rsup|\<ast\>>||2>\<leqslant\><frac|\<lambda\><rsub|j>|\<omega\>>(<enorm|v<rsub|j>||2>-<enorm|K<rsub|j><rsup|*\<ast\>>||2><enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>)>.
    From (<reference|eq:vcycle-smoothingprop>), we conclude

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|K<rsub|j>v<rsub|j>||2>>|<cell|\<leqslant\>>|<cell|a((id-\<theta\><rsub|j>A<rsub|j>)v<rsub|j>,v<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|<enorm|v<rsub|j>||2>-<frac|\<omega\>|\<lambda\><rsub|j>><enorm|A<rsub|j>v<rsub|j>||2>,>>>>
    </eqnarray*>

    which implies, letting <with|mode|math|v<rsub|j>\<rightarrow\>K<rsub|j><rsup|\<ast\>>v<rsub|j>>
    and <with|mode|math|<wide|K|\<bar\>>\<assign\>K<rsub|j>K<rsub|j><rsup|\<ast\>>>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<frac|\<omega\>|\<lambda\><rsub|j>><enorm|A<rsub|j>K<rsup|\<ast\>><rsub|j>v<rsub|j>||2>>|<cell|\<leqslant\>>|<cell|<enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>-<enorm|K<rsub|j>K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>>>|<row|<cell|>|<cell|=>|<cell|a(K<rsub|j><rsup|\<ast\>>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)-a(<wide|K|\<bar\>>v<rsub|j>,<wide|K|\<bar\>>v<rsub|j>)=a((<wide|K|\<bar\>><rsub|j>-<wide|K|\<bar\>><rsub|j><rsup|2>)v<rsub|j>,v<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|a(<wide|K|\<bar\>>(id-<wide|K|\<bar\>>)v<rsub|j>,v<rsub|j>)\<leqslant\>a((id-<wide|K|\<bar\>>)v<rsub|j>,v<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|<enorm|v<rsub|j>||2>-<enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>.>>>>
    </eqnarray*>

    <em|Step 4.> Show

    <\equation*>
      a(E<rsub|j>v<rsub|j>,v<rsub|j>)\<leqslant\>\<delta\><enorm|v<rsub|j>||2>
    </equation*>

    by induction. For <with|mode|math|j=0>, <with|mode|math|E<rsub|j>=0>:
    easy. For <with|mode|math|j-1\<rightarrow\>j>:

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(E<rsub|j>v<rsub|j>,v<rsub|j>)>|<cell|=>|<cell|a(K<rsub|j>(id-P<rsub|j-1>+E<rsub|j-1>P<rsub|j-1>)K<rsup|\<ast\>><rsub|j>v<rsub|j>,v<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)+a(E<rsub|j-1>P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>,P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>)>>|<row|<cell|>|<cell|<below|\<leqslant\>|<with|mode|text|Ind.Ass.>>>|<cell|a((id-P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)+\<delta\>a(P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>,P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)+\<delta\>a(P<rsub|j-1>K<rsub|j><rsup|\<ast\>>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(1-\<delta\>)a((id-P<rsub|j-1>)K<rsup|\<ast\>><rsub|j>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)+\<delta\>a(K<rsub|j><rsup|\<ast\>>,v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)>>|<row|<cell|>|<cell|<below|\<leqslant\>|(\<ast\>)>>|<cell|(1-\<delta\>)<frac|C<rsub|p>|\<omega\>>(<enorm|v<rsub|j>||2>-<enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>)+\<delta\><enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j><rsup|\<ast\>>||2>=\<delta\><enorm|v<rsub|j>||2>.>>>>
    </eqnarray*>

    To prove <with|mode|math|(\<ast\>)>, consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|a((id-P<rsub|j-1>)K<rsup|\<ast\>><rsub|j>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)<rsup|2>>|<cell|=>|<cell|<ip|(id-P<rsub|j-1>)K<rsup|\<ast\>><rsub|j>v<rsub|j>|A<rsub|j>K<rsub|j><rsup|\<ast\>>v<rsub|j>||><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|(id-P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>v<rsub|j>||2><l2norm|A<rsub|j>K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>>>|<row|<cell|>|<cell|=>|<cell|<l2norm|(id-P<rsub|j-1>)<rsup|2>K<rsub|j><rsup|\<ast\>>v<rsub|j>||2><l2norm|A<rsub|j>K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>>>|<row|<cell|>|<cell|=>|<cell|C<rsub|p>\<lambda\><rsub|j><rsup|-1><wide*|<enorm|(id-P<rsub|j-1>)K<rsub|j><rsup|\<ast\>>v<rsub|j>||2>|\<wide-underbrace\>><rsub|=a((id-P<rsub|j-1>)K<rsup|\<ast\>><rsub|j>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)>,>>>>
    </eqnarray*>

    which implies

    <\eqnarray*>
      <tformat|<table|<row|<cell|a((id-P<rsub|j-1>)K<rsup|*\<ast\>><rsub|j>v<rsub|j>,K<rsub|j><rsup|\<ast\>>v<rsub|j>)>|<cell|\<leqslant\>>|<cell|C<rsub|p><frac|1|\<omega\>><left|(><enorm|v<rsub|j>||2>-<enorm|K<rsub|j><rsup|\<ast\>>v<rsub|j>||2><right|)>.>>>>
    </eqnarray*>

    Also,

    <\equation*>
      (1-\<delta\>)<frac|C<rsub|p>|\<omega\>>=<left|(>1-<frac|C<rsub|p>|C<rsub|p>+\<omega\>><right|)><frac|C<rsub|p>|\<omega\>>=<frac|\<omega\>|C<rsub|p>+\<omega\>>*<frac|C<rsub|p>|\<omega\>>=<frac|C<rsub|p>|C<rsub|p>+w>=\<delta\>.
    </equation*>

    \;
  </proof>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|5|?>>
    <associate|auto-2|<tuple|5.1|?>>
    <associate|eq:smoothingprop|<tuple|5.1|?>>
    <associate|eq:vcycle-approxprop|<tuple|5.2|?>>
    <associate|eq:vcycle-smoothingprop|<tuple|5.1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>The
      symmetric V-cycle> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>