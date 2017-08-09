export function setupXLabel(svg, width, height, label, xAxis) {
  svg.select('.x.axis')
    .attr('transform', `translate(0,${height})`)
    .transition()
    .call(xAxis);

  svg.select('.x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .text(label)
    .attr('transform', `translate(${(width / 2)},${(height + 30)})`);
}

export function setupYLabel(svg, height, label, yAxis) {
  // for some reason using transition here breaks the text calculation below
  const a = svg.select('.y.axis')
    .call(yAxis);

  const l = svg.select('.y.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .text(label);
  const all = a.selectAll('text')._groups[0];
  const textHeight = all[all.length - 1].getComputedTextLength() + 20;
  l.attr('transform', `translate(-${textHeight},${(height / 2)})rotate(-90)`);
  return textHeight + 20;
}
