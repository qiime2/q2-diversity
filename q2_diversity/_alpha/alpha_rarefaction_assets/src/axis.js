export function setupXLabel(svg, width, height, label, xAxis) {
  let maxLabelX = 30;
  svg.select('.x.axis')
    .attr('transform', `translate(0,${height})`)
    .call(xAxis)
      .selectAll('text')
      .style('text-anchor', 'end')
      .attr('dx', '-.8em')
      .attr('dy', '-0.5em')
      .attr('transform', function setHeight() {
        const textHeight = this.getComputedTextLength();
        if (textHeight > maxLabelX) maxLabelX = textHeight;
        return 'rotate(-90)';
      });

  svg.select('.x.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(500,${(maxLabelX + height)})`)
    .text(label);
  return maxLabelX;
}

export function setupYLabel(svg, height, label, yAxis) {
  let maxLabelY = 10;
  svg.select('.y.axis')
    .call(yAxis)
      .selectAll('text')
      .attr('transform', function setHeight() {
        const textHeight = this.getComputedTextLength();
        if (textHeight > maxLabelY) maxLabelY = textHeight;
        return '';
      });

  svg.select('.y.label')
    .attr('text-anchor', 'middle')
    .style('font', '12px sans-serif')
    .attr('transform', `translate(-${maxLabelY},${(height / 2)})rotate(-90)`)
    .text(label);
  return maxLabelY;
}
