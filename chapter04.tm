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

  \;
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
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>The
      W-Cycle> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>