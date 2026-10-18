"""Display-only helpers; preserve correlations, tests and MDS coordinates."""
import numpy as np
from matplotlib.transforms import Bbox

def ranking_text(rank,tests):
    rows=list(rank.sort_values('rho_with_tD1',ascending=False).itertuples())
    parts=[]
    for i,r in enumerate(rows):
        if i:
            prev=rows[i-1].comparator
            t=tests[((tests.A==prev)&(tests.B==r.comparator))|((tests.B==prev)&(tests.A==r.comparator))]
            assert len(t)==1
            parts.append(' ≈ (n.s.) ' if t.iloc[0].p_unadjusted>=.05 else ' > ')
        parts.append(f'{r.comparator} ({r.rho_with_tD1:.4f})')
    return ''.join(parts)

def place_labels(fig,ax,xy,names,colors,fontsize=14):
    """Choose labels inside axes with no text-text or text-marker box overlap."""
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    points=ax.transData.transform(xy)
    bounds=ax.get_window_extent().padded(-7)
    obstacles=[Bbox.from_bounds(x-10,y-10,20,20) for x,y in points]
    texts=[ax.text(0,0,n,ha='center',va='center',fontsize=fontsize,color=c,zorder=5) for n,c in zip(names,colors)]
    fig.canvas.draw()
    sizes=[(t.get_window_extent(renderer).width+10,t.get_window_extent(renderer).height+8) for t in texts]
    candidates=[]
    for (px,py),(w,h) in zip(points,sizes):
        opts=[]
        for gap in [1,4,8,14,22,34]:
            for dx,dy in [(0,1),(0,-1),(1,0),(-1,0),(1,1),(-1,1),(1,-1),(-1,-1)]:
                cx=px+dx*(w/2+10+gap);cy=py+dy*(h/2+10+gap)
                b=Bbox.from_bounds(cx-w/2,cy-h/2,w,h)
                if bounds.contains(b.x0,b.y0) and bounds.contains(b.x1,b.y1) and not any(b.overlaps(o) for o in obstacles):opts.append((cx,cy,b))
        candidates.append(opts)
    order=sorted(range(len(names)),key=lambda i:len(candidates[i]));chosen={}
    def solve(k):
        if k==len(order):return True
        i=order[k]
        for candidate in candidates[i]:
            if not any(candidate[2].overlaps(v[2]) for v in chosen.values()):
                chosen[i]=candidate
                if solve(k+1):return True
                del chosen[i]
        return False
    if not solve(0):raise RuntimeError('No collision-free label layout')
    for i,t in enumerate(texts):
        cx,cy,b=chosen[i];t.set_position(ax.transData.inverted().transform([cx,cy]))
    fig.canvas.draw()
    boxes=[t.get_window_extent(fig.canvas.get_renderer()).padded(3) for t in texts]
    assert all(bounds.contains(b.x0,b.y0) and bounds.contains(b.x1,b.y1) for b in boxes)
    assert all(not b.overlaps(o) for b in boxes for o in obstacles)
    assert all(not boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i))
