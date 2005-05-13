<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|1><section|A Two-Level Method>

  <subsection|The Model Problem>

  We consider a weak formulation of the Laplace problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|-\<Delta\>u>|<cell|=>|<cell|f<space|1fn><with|mode|text|in
    <with|mode|math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>>>,>>|<row|<cell|u>|<cell|=>|<cell|0<space|1fn><with|mode|text|on
    <with|mode|math|\<partial\>\<Omega\>>>.>>>>
  </eqnarray*>

  Let <with|mode|math|v:\<Omega\>\<rightarrow\>\<bbb-R\>> be a test function
  with <with|mode|math|v\|<rsub|\<partial\>\<Omega\>>=0>. Gauÿ:
  <with|mode|math|\<sigma\>\<assign\>v\<nabla\>u>

  <\equation*>
    <big|int><rsub|\<Omega\>>div \<sigma\>*d
    v=<big|int><rsub|\<partial\>\<Omega\>>\<sigma\>\<cdot\>n*d a=0
  </equation*>

  <\equation*>
    div \<sigma\>=\<nabla\>v\<cdot\>\<nabla\>u+v\<Delta\>u=\<nabla\>u\<cdot\>\<nabla\>v-f*v
  </equation*>

  <\equation*>
    \<Rightarrow\><big|int><rsub|\<Omega\>>\<nabla\>u\<cdot\>\<nabla\>v=<big|int>f*v*d
    v.
  </equation*>

  <\definition>
    Let

    <\eqnarray*>
      <tformat|<table|<row|<cell|<ip|v|w>>|<cell|\<assign\>>|<cell|<big|int><rsub|\<Omega\>>v*w*d
      v<space|1fn><with|mode|text|(inner product in
      <with|mode|math|L<rsup|2>(\<Omega\>)>)>,>>|<row|<cell|a(v,w)>|<cell|\<assign\>>|<cell|<big|int><rsub|\<Omega\>>\<nabla\>v\<cdot\>\<nabla\>w*d
      v<space|1fn><with|mode|text|(<em|energy inner product> in
      <with|mode|math|L<rsup|2>(\<Omega\>)>)>,>>>>
    </eqnarray*>

    and <with|mode|math|\<\|\|\>v\<\|\|\>\<assign\><sqrt|<ip|v|v>>>,
    <with|mode|math|\<interleave\>v\<interleave\>\<assign\><sqrt|a(v,v)>>.

    <\equation*>
      H<rsup|1>(\<Omega\>)\<assign\>{v:v\<in\>L<rsup|2>(\<Omega\>),\<partial\><rsub|1>v,\<partial\><rsub|2>v\<in\>L<rsup|2>(\<Omega\>)}.
    </equation*>
  </definition>

  <\theorem>
    <with|mode|math|H<rsup|1>(\<Omega\>)> is a Hilbert space with
    <with|mode|math|\<\|\|\>u\<\|\|\><rsub|<rsub|1>>\<assign\><sqrt|\<\|\|\>u\<\|\|\><rsup|2>+\<interleave\>u\<interleave\><rsup|2>>>.
    <with|mode|math|H<rsup|1><rsub|0>(\<Omega\>)\<assign\>clos(C<rsub|0><rsup|\<infty\>>(\<Omega\>),<l2norm|\<cdot\>><rsub|1>)>
    is a Hilbert space with <with|mode|math|<enorm|\<cdot\>>>.
  </theorem>

  <\definition>
    Let <with|mode|math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>> be a polygonal
    domain, <with|mode|math|h\<gtr\>0> the mesh size parameter. A
    <em|uniform, consistent triangulation> <with|mode|math|\<cal-C\><rsub|h>>
    is a decomposition on <with|mode|math|<wide|\<Omega\>|\<bar\>>=<big|cup><rsub|C\<in\>\<cal-C\><rsub|h>><wide|\<Omega\>|\<bar\>><rsub|C>>
    such that

    <big-figure|<postscript|tri-transform.fig|*5/8|*5/8||||>|The triangle
    transformation.>

    <\enumerate-alpha>
      <item><with|mode|math|<wide|\<Omega\>|\<bar\>><rsub|C>=T<rsub|C>(<wide|\<Omega\>|^>)>
      with <with|mode|math|<wide|\<Omega\>|^>=conv(<wide|z<rsub|0>|^>,<wide|z<rsub|1>|^>,<wide|z<rsub|2>|^>)>
      with

      <\eqnarray*>
        <tformat|<table|<row|<cell|T<rsub|C>(<wide|z|^>)>|<cell|\<assign\>>|<cell|z<rsub|0>+J<rsub|C><wide|z|^>,>>|<row|<cell|D*T<rsub|C>>|<cell|=>|<cell|J<rsub|C>=<matrix|<tformat|<cwith|1|1|1|1|cell-halign|r>|<table|<row|<cell|z<rsub|1>-z<rsub|0>>|<cell|z<rsub|2>-z<rsub|0>>>>>>,>>>>
      </eqnarray*>

      <item><with|mode|math|det J<rsub|C>\<gtr\>0>,
      <with|mode|math|\|J<rsub|C>\|\<leqslant\>C<rsub|\<Omega\>>h>,
      <with|mode|math|\|J<rsub|C><rsup|-1>\|\<leqslant\>C<rsub|\<Omega\>>h<rsup|-1>>
      with <with|mode|math|\|J\|\<assign\>sup<rsub|\|<wide|z|^>\|=1>\|J*<wide|z|^>\|>,
      <with|mode|math|\|<wide|z|^>\|\<assign\><sqrt|<wide|z|^><rsub|1><rsup|2>+<wide|z<rsub|2>|^><rsup|2>>>,
      <with|mode|math|z<rsub|i>=T<rsub|C><wide|z|^><rsub|i>> (mesh
      quasi-uniformness),

      <item><with|mode|math|<wide|\<Omega\>|\<bar\>><rsub|C>\<cap\><wide|\<Omega\>|\<bar\>><rsub|C<rprime|'>>=conv({z<rsub|0>,z<rsub|1>,z<rsub|2>}\<cap\>{z<rsub|0><rprime|'>,z<rsub|1><rprime|'>,z<rsub|2><rprime|'>}>
      (mesh admissibility).
    </enumerate-alpha>
  </definition>

  <\theorem>
    <\enumerate-alpha>
      <item><with|mode|math|X<rsub|h>\<assign\>{v\<in\>C(<wide|\<Omega\>|\<bar\>>):v<rsub|h>\|<rsub|<wide|\<Omega\>|\<bar\>><rsub|C>><with|mode|text|
      linear>}\<subset\>H<rsup|1>(\<Omega\>)>.

      <item><with|mode|math|v<rsub|h>\<in\>X<rsub|h>> is uniquely defined by
      the nodal values <with|mode|math|v<rsub|h>(z)> at
      <with|mode|math|z\<in\><wide|\<cal-N\>|\<bar\>><rsub|h>\<assign\><big|cup><rsub|C\<in\>\<cal-C\><rsub|h>><big|cup><rsub|i>T<rsub|C>(<wide|z|^><rsub|i>)>.

      <item><with|mode|math|V<rsub|h>\<assign\>X<rsub|h>\<cap\>H<rsup|1><rsub|0>(\<Omega\>)=span{\<varphi\><rsub|z>\<in\>X<rsub|h>:z\<in\>\<cal-N\><rsub|h>\<assign\><wide|\<cal-N\>|\<bar\>><rsub|h>\<setminus\>\<partial\>\<Omega\>}>
      with

      <\equation*>
        \<varphi\><rsub|z>(y)=<choice|<tformat|<cwith|1|1|1|1|cell-valign|b>|<table|<row|<cell|1>|<cell|y=z>>|<row|<cell|0>|<cell|<with|mode|text|elsewhere>>>>>><space|1fn>(y\<in\>\<cal-N\><rsub|h>).
      </equation*>
    </enumerate-alpha>
  </theorem>

  <\lemma>
    The finite element problem of finding
    <with|mode|math|u<rsub|h>\<in\>V<rsub|h>:><with|mode|math|a(u<rsub|h>,v<rsub|h>)=(f,v<rsub|h>)>
    for all <with|mode|math|v\<in\>V<rsub|h>> has a unique solution.
  </lemma>

  <\proof>
    <with|mode|math|\<b-A\>=(a(\<varphi\><rsub|z>,\<varphi\><rsub|y>))<rsub|z,y\<in\>\<cal-N\><rsub|h>>>,
    <with|mode|math|\<b-f\>=((f,y<rsub|z>))<rsub|z\<in\>\<cal-N\><rsub|h>>>.
    <with|mode|math|\<b-u\>=\<b-A\><rsup|-1>\<b-f\>>, where
    <with|mode|math|u<rsub|h>=<big|sum><rsub|z\<in\>\<cal-N\><rsub|h>>\<b-u\>[z]\<varphi\><rsub|z>>.
    Lax-Milgram ensures the existence of <with|mode|math|\<b-A\><rsup|-1>>.
  </proof>

  <\example*>
    <with|mode|math|\<Omega\>=(0,1)<rsup|2>>, <with|mode|math|h=1/n>. Then
    <with|mode|math|\<b-A\>> is exactly the same as in Section
    <reference|subsec:fd-model-problem>.

    We obtain <with|mode|math|\|J<rsub|C>\|=h>,
    <with|mode|math|\|J<rsub|C><rsup|-1>\|=h<rsup|-1>><with|mode|math|\<Rightarrow\>><with|mode|math|C<rsub|\<Omega\>>=1>.

    <big-figure|<postscript|uniform-grid-unit-square.fig|4cm|||||>|A uniform
    grid on <with|mode|math|(0,1)<rsup|2>>.>
  </example*>

  We define <with|mode|math|A<rsub|h>:V<rsub|h>\<rightarrow\>V<rsub|h>> to be
  <with|mode|math|(A<rsub|h>v<rsub|h>,w<rsub|h>)=a(v<rsub|h>,w<rsub|h>)> for
  <with|mode|math|v<rsub|h>,w<rsub|h>\<in\>V<rsub|h>>.

  <\equation*>
    \<b-u\><rsup|k+1>=\<b-u\><rsup|k>+\<theta\>(\<b-f\>-\<b-A\>\<b-u\><rsup|k>)
  </equation*>

  corresponds to

  <\equation*>
    u<rsup|k+1><rsub|h>=u<rsup|k><rsub|h>\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k>)<rsup|h>.
  </equation*>

  The former converges if <with|mode|math|\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=\<\|\|\>\<b-I\>-\<theta\>\<b-A\>\<\|\|\>\<less\>1>.
  (For the equality, remember that <with|mode|math|\<b-A\>> is symmetric.)
  The latter converges if <with|mode|math|\<\|\|\>id-\<theta\><rsub|h>A<rsub|h>\<\|\|\>\<less\>1\<Leftrightarrow\>\<theta\><rsub|h>\<\|\|\>A<rsub|h>\<\|\|\>\<less\>1>,
  where <with|mode|math|\<\|\|\>A<rsub|h>\<\|\|\>\<assign\>sup<rsub|\<\|\|\>v<rsub|h>\<\|\|\>=1>\<\|\|\>A<rsub|h>v<rsub|h>\<\|\|\>>.

  <big-figure|<postscript|red-refinement.fig|*5/8|*5/8||||>|Red refinement.>

  Let <with|mode|math|\<cal-C\><rsub|H>> be a uniform, consistent
  triangulation. Then, define <with|mode|math|\<cal-C\><rsub|h>> with
  <with|mode|math|h=1/2H> by uniform refinement, as in

  <\equation*>
    <wide|\<Omega\>|\<bar\>>=<big|cup><rsub|C\<in\>\<cal-C\><rsub|H>><big|cup><rsub|i=0><rsup|3>T<rsub|C>(<wide|\<Omega\>|^><rsub|i>)
  </equation*>

  with

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|\<Omega\>|^><rsub|0>>|<cell|=>|<cell|conv{<wide|z|^><rsub|0>,<wide|z|^><rsub|01>,<wide|z|^><rsub|02>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|1>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|01>,<wide|z|^><rsub|12>,<wide|z|^><rsub|02>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|2>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|01>,<wide|z|^><rsub|1>,<wide|z|^><rsub|12>},>>|<row|<cell|<wide|\<Omega\>|^><rsub|3>>|<cell|=>|<cell|conv
    {<wide|z|^><rsub|02>,<wide|z|^><rsub|12>,<wide|z|^><rsub|2>},>>>>
  </eqnarray*>

  which implies <with|mode|math|V<rsub|H>\<subset\>V<rsub|h>><emdash>the
  finite element spaces are <em|nested>.

  <\equation*>
    A<rsub|h>u<rsub|h>=f<rsub|h>\<Leftrightarrow\><ip|A<rsub|h>u<rsub|h>|v<rsub|h>>=<ip|f<rsub|h>|v<rsub|h>><space|1fn><with|mode|text|for
    all <with|mode|math|v<rsub|h>\<in\>V<rsub|h>>>.
  </equation*>

  Now define <with|mode|math|c<rsub|H>\<in\>V<rsub|H>> as a form of
  coarse-space residual by requiring

  <\equation*>
    a(c<rsub|H>,v<rsub|H>)=<ip|A<rsub|H>c<rsub|H>|v<rsub|H>><above|=|!><ip|f<rsub|h>-A<rsub|h>u<rsup|k><rsub|h>|v<rsub|H>>=<ip|f<rsub|h>|v<rsub|H>>-a(u<rsup|k>,v<rsub|H>)<space|1fn><with|mode|text|for
    all <with|mode|math|v<rsub|H>\<in\>V<rsub|H>.>>
  </equation*>

  (<with|mode|math|c<rsub|H>> may also be viewed as an element of
  <with|mode|math|V<rsub|h>> through interpolation.)

  <\note>
    <em|(Two-level method)> Let <with|mode|math|u<rsub|h><rsup|0>\<in\>V<rsub|h>>.
    For <with|mode|math|k=1,2,\<ldots\>>

    <\enumerate>
      <item>compute <with|mode|math|u<rsub|h><rsup|k-1/2>=u<rsup|k-1><rsub|h>+\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>)>,

      <item>compute <with|mode|math|c<rsup|k><rsub|H>\<in\>V<rsub|h>>:
      <with|mode|math|a(c<rsub|H><rsup|k>,v<rsub|H>)=<ip|f<rsub|h>|v<rsub|H>>-a(u<rsub|h><rsup|k-1/2>,v<rsub|H>)>,

      <item>set <with|mode|math|u<rsup|k>=u<rsub|h><rsup|k-1/2>+c<rsub|H><rsup|k>>.
    </enumerate>
  </note>

  <\lemma>
    Let

    <\eqnarray*>
      <tformat|<table|<row|<cell|Q<rsub|H>:L<rsup|2>(\<Omega\>)\<rightarrow\>V<rsub|H><with|mode|text|
      be the <with|mode|math|L<rsup|2>>-projection, i.e.
      >>|<cell|>|<cell|<ip|Q<rsub|H>v|w<rsub|H>>=<ip|v|w<rsub|H>>,>>|<row|<cell|P<rsub|H>:H<rsup|1><rsub|0>(\<Omega\>)\<rightarrow\>V<rsub|H><with|mode|text|
      be the Galerkin projection, i.e. >>|<cell|>|<cell|a(P<rsub|H>v,w<rsub|H>)=a(v,w<rsub|H>)>>>>
    </eqnarray*>

    with <with|mode|math|v\<in\>L<rsup|2>(\<Omega\>)\<supset\>V<rsub|h>>,
    <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>,
    <with|mode|math|v\<in\>H<rsup|1><rsub|0>(\<Omega\>)\<supset\>V<rsub|h>>,
    for all <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>. Then, we have for the
    error propagation of the two-level method

    <\equation*>
      e<rsup|k><rsub|h>=<wide*|u<rsup|k><rsub|h>|\<wide-underbrace\>><rsub|=A<rsup|-1><rsub|h>f<rsub|h>>-u<rsub|h>=(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>.
    </equation*>
  </lemma>

  <\proof>
    \;

    <\enumerate-alpha>
      <item>

      <\eqnarray*>
        <tformat|<table|<row|<cell|e<rsup|k-1/2>>|<cell|=>|<cell|u<rsup|k-1/2><rsub|h>-u<rsub|h>>>|<row|<cell|>|<cell|=>|<cell|u<rsup|k-1>-u<rsub|h>+\<theta\><rsub|h>(A<rsub|h>u<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>)=(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>.>>>>
      </eqnarray*>

      <item>

      <\equation*>
        (A<rsub|H>c<rsub|H><rsup|k>,v<rsub|H>)=a(c<rsub|H><rsup|k>,v<rsub|H>)=a(u<rsub|h>,v<rsub|H>)-a(u<rsub|h><rsup|k-1/2>,v<rsub|H>)=<ip|Q<rsub|H>f<rsub|h>|v<rsub|H>>-<ip|A<rsub|h>u<rsup|k-1/2><rsub|h>|v<rsub|H>>=<ip|Q<rsub|H>f<rsub|h>-Q<rsub|H>A<rsub|h>u<rsup|k-1/2><rsub|h>|v<rsub|H>>
      </equation*>

      Then

      <\equation*>
        A<rsub|H>c<rsub|H><rsup|k>=Q<rsub|H>(<wide*|f<rsub|h>|\<wide-underbrace\>><rsub|=A<rsub|h>u<rsub|h>>-A<rsub|h>u<rsup|k-1/2><rsub|h>)=-Q<rsub|H>A<rsub|h>e<rsup|k-1/2><rsub|h>,
      </equation*>

      so

      <\equation*>
        e<rsup|k><rsub|h>=u<rsup|k><rsub|h>-u<rsub|h>=u<rsup|k-1/2><rsub|h>-u<rsub|h>+c<rsub|H>=e<rsup|k-1/2><rsub|h>+A<rsub|H><rsup|-1>(-Q<rsub|H>A<rsub|h>e<rsup|k-1/2><rsub|h>).
      </equation*>

      <item>For <with|mode|math|w<rsub|H>\<in\>V<rsub|H>>,

      <\equation*>
        <ip|A<rsub|H>P<rsub|H>v<rsub|h>|w<rsub|H>>=a(P<rsub|H>v<rsub|h>,w<rsub|H>)=a(v<rsub|h>,w<rsub|H>)=<ip|A<rsub|h>v<rsub|h>|w<rsub|H<with|color|red|>>>=<ip|Q<rsub|H>A<rsub|h>v<rsub|h>|w<rsub|H>>,
      </equation*>

      so

      <\equation*>
        A<rsub|H>P<rsub|H>=Q<rsub|H>A<rsub|h>\<Rightarrow\>P<rsub|H>=A<rsup|-1><rsub|H>Q<rsub|H>A<rsub|h>.
      </equation*>
    </enumerate-alpha>

    \;
  </proof>

  <\remark>
    We observe that

    <\itemize-dot>
      <item><with|mode|math|Q<rsub|H>> is an orthogonal projection.
      <with|mode|math|Q<rsub|H><rsup|2>=Q<rsub|H>>,
      <with|mode|math|Q<rsub|H>(id-Q<rsub|H>)=0>,
      <with|mode|math|<l2norm|Q<rsub|H>>=1>,
      <with|mode|math|<l2norm|P<rsub|H>>\<gtr\>1>!

      <item><with|mode|math|P<rsub|H>> is an orthogonal projection.
      <with|mode|math|P<rsub|H><rsup|2>=P<rsub|H>>,
      <with|mode|math|P<rsub|H>(id-P<rsub|H>)=0>,
      <with|mode|math|<enorm|P<rsub|H>>=1>.
    </itemize-dot>
  </remark>

  <\theorem>
    <label|the:twolevel>Let

    <\enumerate>
      <item><with|mode|math|<enorm|v<rsub|h>><rsup|2>\<leqslant\>\<theta\><rsub|h><rsup|-1><l2norm|v<rsub|h>><rsup|2>>
      (<with|mode|math|\<Leftrightarrow\>><with|mode|math|\<theta\><rsub|h><l2norm|A<rsub|h>>\<leqslant\>1>),

      <item><with|mode|math|<l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>\<leqslant\>C\<theta\><rsub|h><enorm|v<rsub|h>><rsup|2>>
      (<with|mode|math|\<leqslant\>C<l2norm|v<rsub|h>><rsup|2>><with|mode|math|\<Rightarrow\>><with|mode|math|C\<geqslant\>1>).
    </enumerate>

    Then, we have <with|mode|math|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)><rsup|2>\<leqslant\>1-1/C\<less\>1>,
    i.e. <with|mode|math|<enorm|u<rsup|k><rsub|h>-u<rsub|h>>\<leqslant\>(1-1/C)<rsup|k/2><enorm|u<rsub|h><rsup|0>-u<rsub|h>>>.
  </theorem>

  1) and 2) of Theorem <reference|the:twolevel> imply

  <\equation*>
    <l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>\<leqslant\>C<l2norm|v<rsub|h>><rsup|2>\<Rightarrow\>C\<geqslant\>1.
  </equation*>

  The following two lemmata provide the (seemingly strong) assumptions of
  this theorem.

  <\lemma>
    <label|lem:2step-energy-estimate>

    <\enumerate-alpha>
      <item><with|mode|math|<l2norm|\<nabla\><wide|v|^>><rsub|<wide|\<Omega\>|^>>\<leqslant\><sqrt|48><l2norm|<wide|v|^>><rsub|<wide|\<Omega\>|^>>>
      for linear <with|mode|math|<wide|v|^>(x)=(1-<wide|x|^><rsub|1>-<wide|x|^><rsub|2>)v<rsub|0>+<wide|x|^><rsub|1>v<rsub|1>+<wide|x|^><rsub|2>v<rsub|2>>.

      <item><with|mode|math|<l2norm|\<nabla\>v<rsub|h>>\<leqslant\>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>>>
      for <with|mode|math|v<rsub|h>\<in\>X<rsub|h>> with
      <with|mode|math|C<rsub|I>=<sqrt|48>C<rsub|\<Omega\>>>.
    </enumerate-alpha>
  </lemma>

  <with|mode|math|C<rsub|\<Omega\>>> measures the quality of the
  triangulation. Roughly, <with|mode|math|C<rsub|\<Omega\>>\<thickapprox\>1>
  for ``rectangular'' triangles. If using uniform refinement,
  <with|mode|math|C<rsub|\<Omega\>>> is independent of the mesh size.

  <\proof>
    <\enumerate-alpha>
      <item>Using

      <\equation*>
        \<nabla\><wide|v|^>=<matrix|<tformat|<table|<row|<cell|v<rsub|1>-v<rsub|0>>>|<row|<cell|v<rsub|2>-v<rsub|0>>>>>>,
      </equation*>

      and for <with|mode|math|P\<in\>\<bbb-P\><rsub|2>>

      <\equation*>
        <l2norm|P><rsub|<wide|\<Omega\>|^>>=<frac|1|2>\<cdot\><frac|1|3><left|(>P<left|(><frac|1|2>(<wide|z<rsub|0>|^>+<wide|z<rsub|1>|^>)<right|)>+P<left|(><frac|1|2>(<wide|z<rsub|1>|^>+<wide|z<rsub|2>|^>)<right|)>+P<left|(><frac|1|2>(<wide|z<rsub|0>|^>+<wide|z<rsub|2>|^>)<right|)><right|)>
      </equation*>

      and

      <\equation*>
        <frac|1|2>(a+b)<rsup|2>\<leqslant\>a<rsup|2>+b<rsup|2>,
      </equation*>

      we compute

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|\<nabla\><wide|v|^><rsub|h>><rsub|<wide|\<Omega\>|^>><rsup|2>>|<cell|=>|<cell|<frac|1|2>(v<rsub|1>-v<rsub|0>)<rsup|2>+<frac|1|2>(v<rsub|2>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|2>(v<rsub|1>-v<rsub|2>+v<rsub|2>-v<rsub|0>)<rsup|2>+<frac|1|2>(v<rsub|2>-v<rsub|1>+v<rsub|1>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(v<rsub|1>-v<rsub|2>)<rsup|2>+(v<rsub|2>-v<rsub|0>)<rsup|2>+(v<rsub|2>-v<rsub|1>)<rsup|2>+(v<rsub|1>-v<rsub|0>)<rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|2\<cdot\>4\<cdot\>6\<cdot\><frac|1|6>\<cdot\><frac|1|4><left|(>(v<rsub|1>+v<rsub|2>)<rsup|2>+(v<rsub|2>+v<rsub|0>)<rsup|2>+(v<rsub|1>+v<rsub|0>)<rsup|2><right|)>>>|<row|<cell|>|<cell|=>|<cell|48<l2norm|<wide|v|^><rsub|h>><rsub|<wide|\<Omega\>|^>><rsup|2>.>>>>
      </eqnarray*>

      <item>On the transformed triangle, with
      <with|mode|math|x=T<rsub|C>(<wide|x|^>)> and considering
      <with|mode|math|D*v(x)=\<nabla\>v(x)<rsup|T>>,

      <\eqnarray*>
        <tformat|<table|<row|<cell|<wide|v|^>(<wide|x|^>)>|<cell|=>|<cell|v<rsub|h>(x)>>|<row|<cell|\<Rightarrow\>D<wide|v|^>(<wide|x|^>)>|<cell|=>|<cell|D(v\<circ\>T<rsub|C>(<wide|x|^>))=D*v(x)D*T<rsub|C>>>|<row|<cell|\<Rightarrow\>\<nabla\><wide|v|^>(x)>|<cell|=>|<cell|D<wide|v|^>(<wide|x|^>)<rsup|T>=(D*T<rsub|C>)<rsup|T>\<nabla\>v(x).>>>>
      </eqnarray*>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|\<nabla\>v><rsub|\<Omega\><rsub|C>><rsup|2>=<big|int><rsub|\<Omega\><rsub|C>>\|\<nabla\>v\|<rsup|2>d
        x>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|\<nabla\>v\<circ\>T\|<rsup|2>d<wide|x|^>>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|J<rsub|C><rsup|-T>\<nabla\>v\|<rsup|2>d<wide|x|^>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\|J<rsub|C><rsup|-T>\|<rsup|2>\|det
        J<rsub|C>\|<l2norm|\<nabla\><wide|v|^>><rsup|2><rsub|<wide|\<Omega\>|^>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|\<Omega\>>h<rsup|-2>48<l2norm|<wide|v|^>><rsup|2><rsub|\<Omega\><rsub|C>>>>>>
      </eqnarray*>

      So,

      <\equation*>
        <enorm|v<rsub|h>><rsup|2>=<big|sum><rsub|C\<in\>\<cal-C\><rsub|h>><l2norm|\<nabla\>v<rsub|h>><rsub|\<Omega\><rsub|C>>\<leqslant\>48C<rsub|\<Omega\>><rsup|2>h<rsup|-2><big|sum><rsub|C\<in\>\<cal-C\><rsub|h>><l2norm|v<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>.
      </equation*>
    </enumerate-alpha>

    \;
  </proof>

  <\corollary>
    Let <with|mode|math|A<rsub|h>:V<rsub|h>\<rightarrow\>V<rsub|h>> be
    defined by <with|mode|math|<ip|A<rsub|h>v<rsub|h>|w<rsub|h>>=a(v<rsub|h>,w<rsub|h>)>.
    Then

    <\equation*>
      <l2norm|A<rsub|h>>=sup<rsub|<l2norm|v<rsub|h>>=1><l2norm|A<rsub|h>v<rsub|h>>\<leqslant\>C<rsub|I>h<rsup|-2>.
    </equation*>

    \;
  </corollary>

  <\proof>
    \;

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|A<rsub|h>v<rsub|h>>>|<cell|=>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1><ip|A<rsub|h>v<rsub|h>|w<rsub|h>>=sup<rsub|<l2norm|w<rsub|h>>=1>a(v<rsub|h>,w<rsub|h>)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1><enorm|v<rsub|h>><enorm|w<rsub|h>>>>|<row|<cell|>|<cell|<below|\<leqslant\>|(<with|mode|text|<reference|lem:2step-energy-estimate>)>>>|<cell|sup<rsub|<l2norm|w<rsub|h>>=1>C<rsub|I>h<rsup|-1><l2norm|v<rsub|h>>C<rsub|I>h<rsup|-1><l2norm|w<rsub|h>>=C<rsub|I>h<rsup|-2><l2norm|v<rsub|h>>.>>>>
    </eqnarray*>
  </proof>

  <\lemma>
    <\enumerate-alpha>
      <item>Let <with|mode|math|<wide|v|^>\<in\>C(<wide|\<Omega\>|^>)> be
      linear on <with|mode|math|<wide|\<Omega\>|^><rsub|i>>, let
      <with|mode|math|<wide|\<Pi\>|^><wide|v|^>> be linear on
      <with|mode|math|<wide|\<Omega\>|^>> with
      <with|mode|math|<wide|\<Pi\>|^><wide|v|^>(<wide|z|^><rsub|i>)=<wide|v|^>(<wide|z|^><rsub|i>)>.
      Then, we have

      <\equation*>
        <l2norm|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>><rsub|<wide|\<Omega\>|^>>\<leqslant\>3<l2norm|\<nabla\><wide|v|^>><rsub|<wide|\<Omega\>|^>>.
      </equation*>

      <item>

      <\equation*>
        <l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>>\<leqslant\>2<sqrt|3>C<rsub|\<Omega\>>h<enorm|v<rsub|h>>.
      </equation*>
    </enumerate-alpha>
  </lemma>

  <\proof>
    \;

    <\enumerate-alpha>
      <item>Let <with|mode|math|w<rsub|i,j>\<assign\>v<rsub|i,j>-1/2(v<rsub|i>+v<rsub|j>)>.
      Then

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>><rsub|<wide|\<Omega\>|^>>>|<cell|=>|<cell|<frac|1|4>*<frac|1|6><left|(><left|(><frac|w<rsub|0,2>|2><right|)><rsup|2>+<left|(><frac|w<rsub|0,2>+w<rsub|0,1>|2><right|)><rsup|2>+<left|(><frac|w<rsub|0,1>|2><right|)><rsup|2><right|)>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<frac|1|4><big|sum><rsub|i\<less\>j><left|(><frac|w<rsub|i,j>|2><right|)><rsup|2>\<leqslant\><frac|1|4><big|sum><rsub|i\<less\>j><left|[>(v<rsub|i,j>-v<rsub|i>)<rsup|2>+(v<rsub|i,j>-v<rsub|j>)<rsup|2><right|]>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3<l2norm|\<nabla\><wide|v|^>><rsup|2>.>>>>
      </eqnarray*>

      <item>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-\<Pi\><rsub|H>v<rsub|H>><rsub|\<Omega\><rsub|C>><rsup|2>>|<cell|=>|<cell|<big|int><rsub|\<Omega\><rsub|C>>\|v<rsub|h>-\<Pi\><rsub|H>v<rsub|H>\|<rsup|2>d
        x>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|\|<wide|v|^>-<wide|\<Pi\>|^><wide|v|^>\|d*
        <wide|x|^>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3<big|int><rsub|<wide|\<Omega\>|^>>\|det*J<rsub|C>\|<below|\|\<nabla\><wide|v|^>\|<rsup|2><rsub|>|=\|J<rsub|C><rsup|T>\<nabla\>v<rsub|h>\|>d
        x>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3C<rsub|\<Omega\>><rsup|2>h<rsup|2><l2norm|\<nabla\>v<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>.>>>>
      </eqnarray*>

      <\eqnarray*>
        <tformat|<table|<row|<cell|<l2norm|v<rsub|h>-Q<rsub|H>v<rsub|h>><rsup|2>>|<cell|=>|<cell|min<rsub|V<rsub|H>\<in\>V<rsub|H>><l2norm|v<rsub|h>-V<rsub|H>><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|v<rsub|h>-\<Pi\><rsub|H>v<rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|C\<in\>\<cal-C\><rsub|H>><l2norm|v<rsub|h>-\<Pi\>V<rsub|h>><rsup|2><rsub|\<Omega\><rsub|C>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|3C<rsub|\<Omega\>><rsup|2>h<rsup|2><enorm|v<rsub|h>>.>>>>
      </eqnarray*>
    </enumerate-alpha>
  </proof>

  So, <with|mode|math|\<theta\><rsub|h>=h<rsup|2>/C<rsup|2><rsub|I>> and
  <with|mode|math|C*h<rsup|2>/C<rsup|2><rsub|I>=C<rsup|2><rsub|Q>h<rsup|2>><with|mode|math|\<Rightarrow\>><with|mode|math|C=C<rsup|2><rsub|Q>C<rsup|2><rsub|I>=48\<cdot\>12\<cdot\>C<rsup|4><rsub|\<Omega\>>>.
  We obtain <with|mode|math|\<theta\><rsub|h>=h<rsup|2>/<sqrt|48C<rsub|\<Omega\>>>>.
  <with|color|red|Was um alles in der Welt ist <with|mode|math|C<rsub|Q>>?>
  (s. Wieners 12 unten)

  <\proof>
    (of <reference|the:twolevel>) <with|mode|math|e<rsub|h><rsup|k>=(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>>,\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<rho\>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)>=sup<rsub|<l2norm|v<rsub|h>>><enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)>>>|<row|<cell|\<Leftrightarrow\><enorm|e<rsup|k><rsub|h>>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)<wide*|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>|\<wide-underbrace\>><rsub|e<rsup|k-1/2><rsub|h>>>\<leqslant\>\<rho\><enorm|e<rsup|k-1><rsub|h>>>>>>
    </eqnarray*>

    for all <with|mode|math|e<rsup|k-1><rsub|h>\<in\>V<rsub|h>>.

    <em|1st step.> For all <with|mode|math|v<rsub|H>\<in\>V<rsub|H>>,

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(e<rsup|k><rsub|h>,v<rsub|H>)>|<cell|=>|<cell|a((id-P<rsub|H>)e<rsup|k-1/2><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)<rsup|2>e<rsup|k-1/2><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)e<rsup|k><rsub|h>,v<rsub|H>)>>|<row|<cell|>|<cell|=>|<cell|a(e<rsup|k><rsub|h>,(id-P<rsub|H>)v<rsub|H>)=a(e<rsup|k><rsub|h>,v<rsub|H>)-a(e<rsup|k><rsub|h>,<wide*|P<rsub|H>v<rsub|H>|\<wide-underbrace\>><rsub|=v<rsub|H>>)=0.>>>>
    </eqnarray*>

    <em|2nd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>>|<cell|=>|<cell|a(e<rsup|k><rsub|h>,e<rsup|k><rsub|h>)=a(e<rsup|k><rsub|h>,e<rsup|k><rsub|h>-Q<rsub|H>e<rsup|k><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|A<rsub|h>e<rsup|k><rsub|h>|e<rsub|h><rsup|k>-Q<rsub|H>e<rsup|k><rsub|h>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<l2norm|A<rsub|h>e<rsup|k><rsub|h>><l2norm|e<rsub|h><rsup|k>-Q<rsub|H>e<rsup|k><rsub|h>>\<leqslant\><l2norm|A<rsub|h>e<rsub|h><rsup|k>><sqrt|C\<theta\><rsub|h>><enorm|e<rsup|k><rsub|h>>,>>>>
    </eqnarray*>

    using assumption 2.

    <\equation>
      <label|eq:twolevel-step2>\<Rightarrow\><enorm|e<rsup|k><rsub|h>><rsup|2>\<leqslant\><sqrt|C\<theta\><rsub|h>><l2norm|A<rsub|h>e<rsub|h><rsup|k>>\<Rightarrow\><enorm|e<rsub|h><rsup|k>><rsup|2>\<leqslant\>C\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>
    </equation>

    \;

    <em|3rd step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<Rightarrow\><enorm|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>><rsup|2>>|<cell|=>|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>-2a(e<rsup|k><rsub|h>,\<theta\><rsub|h>A<rsub|h>e<rsup|k><rsub|h>)+<enorm|\<theta\><rsub|h>A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<enorm|e<rsub|h><rsup|k>><rsup|2>-2\<theta\><rsub|h><l2norm|A<rsub|h>e<rsub|h><rsup|k>><rsup|2>+\<theta\><rsub|h><rsup|2><enorm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><l2norm|A<rsub|h>e<rsub|h><rsup|k>><rsup|2>-\<theta\><rsub|h><left|(><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><enorm|A<rsub|h>e<rsup|h><rsub|k>><rsup|2><right|)>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|<below|\<leqslant\>|(<reference|eq:twolevel-step2>)>>|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>-<frac|1|C><enorm|e<rsup|k><rsub|h>><rsup|2>=<left|(>1-<frac|1|C><right|)><enorm|e<rsup|k><rsub|h>><rsup|2>,>>>>
    </eqnarray*>

    where, considering (<reference|eq:twolevel-step2>), we used

    <\equation*>
      -\<theta\><rsub|h><l2norm|A<rsub|h>e<rsup|k><rsub|h>><rsup|2>\<leqslant\>-<frac|1|C><enorm|e<rsup|k><rsub|h>>.
    </equation*>

    <em|4th step.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|e<rsup|k><rsub|h>><rsup|2>>|<cell|=>|<cell|<enorm|(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|a((id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>,(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|a((id-\<theta\><rsub|h>A<rsub|h>)(id-P<rsub|H>)(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k-1><rsub|h>,e<rsup|k-1><rsub|h>)>>|<row|<cell|>|<cell|=>|<cell|<ip|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>|e<rsup|k><rsub|h>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<enorm|(id-\<theta\><rsub|h>A<rsub|h>)e<rsup|k><rsub|h>><enorm|e<rsup|k-1><rsub|h>>\<leqslant\><rsub|3rd
      step><sqrt|1-1/C><enorm|e<rsup|k><rsub|h>><enorm|e<rsup|k-1><rsub|h>>.>>>>
    </eqnarray*>

    \;
  </proof>

  <subsection|Implementation of the Two-Level Method>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A<rsub|H>>|<cell|:>|<cell|V<rsub|H>\<rightarrow\>V<rsub|H>,>>|<row|<cell|>|<cell|>|<cell|(A<rsub|H>v<rsub|H>,w<rsub|H>)=a(v<rsub|H>,w<rsub|H>),>>|<row|<cell|\<b-A\><rsub|H>>|<cell|\<assign\>>|<cell|(a(\<varphi\><rsub|z><rsup|H>,\<varphi\><rsub|y><rsup|H>))<rsub|z,y\<in\>\<cal-N\><rsub|H>>\<in\>\<bbb-R\><rsup|N<rsub|H>\<times\>N<rsub|H>><with|mode|text|
    with >N<rsub|H>:=\|\<cal-N\><rsub|H>\|,>>|<row|<cell|>|<cell|>|<cell|\<b-v\><rsub|H><rsup|T>\<b-A\><rsub|H>\<b-w\><rsub|H>=a(v<rsub|H>,w<rsub|H>),>>|<row|<cell|<with|mode|text|for
    >\<b-v\><rsub|H>>|<cell|=>|<cell|<big|sum>\<b-v\><rsub|H>[z]\<varphi\><rsub|z><rsup|H>\<Leftrightarrow\>\<b-v\><rsub|H>=(v<rsub|H>(z))<rsub|z\<in\>\<cal-N\><rsub|H>>,>>|<row|<cell|\<b-A\><rsub|h>>|<cell|\<in\>>|<cell|\<bbb-R\><rsup|N<rsub|h>\<times\>N<rsub|h>>,>>|<row|<cell|\<b-v\><rsub|H>>|<cell|\<in\>>|<cell|V<rsub|H>\<subset\>V<rsub|h>,>>|<row|<cell|\<b-v\><rsub|h>>|<cell|\<assign\>>|<cell|(v<rsub|H>(z))<rsub|z\<in\>\<cal-N\><rsub|h>>,>>|<row|<cell|>|<cell|>|<cell|\<Rightarrow\>\<b-v\><rsub|n>(z)=<frac|1|2><left|(>\<b-v\><rsub|h>(x)+\<b-v\><rsub|h>(y)<right|)>,>>|<row|<cell|\<b-v\><rsub|h>>|<cell|=>|<cell|\<b-I\><rsub|h><with|mode|text|
    with >\<b-I\><rsub|h>\<in\>\<bbb-R\><rsup|N<rsub|h>\<times\>N<rsub|H>>,>>>>
  </eqnarray*>

  <big-figure|<postscript|2step-refinement.fig|7cm|||||>|Interpolation to a
  finer discretization.>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A<rsub|h>u<rsub|h>=f<rsub|h>>|<cell|\<Leftrightarrow\>>|<cell|<ip|A<rsub|h>u<rsub|h>|\<varphi\><rsub|z><rsup|h>>=<ip|f<rsub|h>|\<varphi\><rsub|z><rsup|h>><space|1fn>\<forall\>z\<in\>\<cal-N\><rsub|h>>>|<row|<cell|>|<cell|\<Leftrightarrow\>>|<cell|\<b-A\><rsub|h>\<b-u\><rsub|h>=\<b-f\><rsub|h><with|mode|text|
    with >\<b-f\><rsub|h>=(<ip|f<rsub|h>|\<varphi\><rsub|h><rsup|z>>)<rsub|z\<in\>\<cal-N\><rsub|h>>>>>>
  </eqnarray*>

  Within finite element context,

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<rsub|h><rsup|k-1/2>>|<cell|=>|<cell|u<rsub|h><rsup|k-1>+\<theta\><rsub|h>(f<rsub|h>-A<rsub|h>u<rsub|h><rsup|k-1>),>>|<row|<cell|<with|mode|text|solve
    >a(c<rsub|H>,\<varphi\><rsub|z><rsup|H>)>|<cell|=>|<cell|<ip|f<rsub|h>|\<varphi\><rsub|z><rsup|H>>-a(u<rsup|k-1/2><rsub|h>,\<varphi\><rsub|z><rsup|H>),>>|<row|<cell|<ip|A<rsub|H>c<rsub|H>|\<varphi\><rsub|z><rsup|H>>>|<cell|=>|<cell|<ip|f<rsub|h>-A<rsub|h>u<rsup|k-1/2><rsub|h>|\<varphi\><rsub|z><rsup|H>>.>>>>
  </eqnarray*>

  Within vector context,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-u\><rsub|h><rsup|k-1/2>>|<cell|=>|<cell|\<theta\><rsub|h>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1>),>>|<row|<cell|(\<b-A\><rsub|H>\<b-c\><rsub|H>)<rsup|T>\<b-v\><rsub|H>>|<cell|=>|<cell|(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>)(\<b-I\><rsub|h>\<b-v\><rsub|H>),>>|<row|<cell|\<Rightarrow\>\<b-v\><rsup|T><rsub|H>\<b-A\><rsub|H>\<b-c\><rsub|H>>|<cell|=>|<cell|\<b-v\><rsub|H><rsup|T>\<b-I\><rsub|h><rsup|T>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>),>>|<row|<cell|\<Rightarrow\>\<b-A\><rsub|H>\<b-c\><rsub|H>>|<cell|=>|<cell|\<b-I\><rsub|h><rsup|T>(\<b-f\><rsub|h>-\<b-A\><rsub|h>\<b-u\><rsub|h><rsup|k-1/2>).>>>>
  </eqnarray*>

  In particular,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<ip|A<rsub|h>v<rsub|H>|w<rsub|H>>>|<cell|=>|<cell|<ip|A<rsub|H>v<rsub|H>|w<rsub|H>>,>>|<row|<cell|\<Leftrightarrow\>(\<b-A\><rsub|h>\<b-I\><rsub|h>\<b-v\><rsub|H>)<rsup|T>(\<b-I\>\<b-w\><rsub|H>)>|<cell|=>|<cell|(\<b-A\><rsub|H>\<b-v\><rsub|H>)<rsup|T>\<b-w\><rsub|H>,>>|<row|<cell|(\<b-v\><rsub|H>)<rsup|T>\<b-I\><rsub|h><rsup|T>>|<cell|=>|<cell|\<b-A\><rsub|h>\<b-I\><rsub|h>\<b-w\><rsub|H>,>>|<row|<cell|\<Rightarrow\>\<b-A\><rsub|H>>|<cell|=>|<cell|\<b-I\><rsub|h><rsup|T>\<b-A\><rsub|h>\<b-I\><rsub|h><with|mode|text|
    (called the <em|Galerkin product>)>.>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|language|german>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|2|?>>
    <associate|auto-2|<tuple|2.1|?>>
    <associate|auto-3|<tuple|2.1|?>>
    <associate|auto-4|<tuple|2.2|?>>
    <associate|auto-5|<tuple|2.3|?>>
    <associate|auto-6|<tuple|2.2|?>>
    <associate|auto-7|<tuple|2.4|?>>
    <associate|eq:twolevel-step2|<tuple|2.1|?>>
    <associate|lem:2step-energy-estimate|<tuple|2.10|?>>
    <associate|the:twolevel|<tuple|2.9|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|The triangle transformation.|<pageref|auto-3>>

      <tuple|normal|A uniform grid on <with|mode|<quote|math>|(0,1)<rsup|2>>.|<pageref|auto-4>>

      <tuple|normal|Red refinement.|<pageref|auto-5>>

      <tuple|normal|Interpolation to a finer
      discretization.|<pageref|auto-7>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|18<space|2spc>A
      Two-Level Method> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|18.1<space|2spc>The Model Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|18.2<space|2spc>Implementation of the
      Two-Level Method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>