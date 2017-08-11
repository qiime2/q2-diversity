import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory20,
  select,
  line,
  nest,
} from 'd3';

import { setupXLabel, setupYLabels } from './axis';
import appendLegendKey from './legend';
import { curData, appendSeries, toggle } from './data';

// re-render chart and legend whenever selection changes
function renderPlot(svg, data, x, y, subY, category, legend, legendTitle) {
  const chart = svg.select('#chart');
  const subChart = svg.select('#subChart');
  const legendBox = select(legend.node().parentNode);

  const depthIndex = data.data.columns.indexOf('depth');
  const medianIndex = data.data.columns.indexOf('50%');
  const countIndex = data.data.columns.indexOf('count');
  let groupIndex = data.data.columns.indexOf('sample-id');
  if (groupIndex === -1) {
    groupIndex = data.data.columns.indexOf(category);
  }
  const points = [data.data.data][0];
  const setGroups = new Set(Array.from(points, d => d[groupIndex]));
  const color = scaleOrdinal(schemeCategory20)
    .domain(setGroups);
  const arrGroups = Array.from(setGroups);

  legend.selectAll('.legend').remove();
  legendTitle.selectAll('.legend').remove();
  legend.attr('height', arrGroups.length * 20);

  let ly = 0;
  const all = 'Select%20All';
  appendSeries(all, [], 'black');
  toggle(all, 'white', null);
  appendLegendKey(legendTitle, all, 10, color);
  const sortedGroupEntries = arrGroups.sort((a, b) => {
    const aIsNaN = isNaN(a);
    const bIsNaN = isNaN(b);
    if (aIsNaN && bIsNaN) {
      // a and b are both alphabetical
      return (a > b) ? 1 : -1;
    } else if (!aIsNaN && bIsNaN) {
      // a is numeric, b is alphabetical
      return 1;
    } else if (aIsNaN && !bIsNaN) {
      // a is alphabetic, b is numeric
      return -1;
    }
    // a and be are both numeric
    return a - b;
  });
  for (const [i, entry] of sortedGroupEntries.entries()) {
    ly = (i + 0.5) * 20;
    const subset = points.filter(d => d[groupIndex] === entry)
                    .sort((a, b) => a[depthIndex] - b[depthIndex]);
    const curColor = color(entry);
    appendSeries(entry, subset, curColor);
    toggle(entry, 'white', null);
    appendLegendKey(legend, entry, ly, color);
  }
  // DOTS
  function plotDots(selection, index, yScale) {
    selection.transition()
      .attr('class', d => `circle ${d[groupIndex]}`)
      .attr('fill', d => color(d[groupIndex]))
      .attr('opacity', d => curData[d[groupIndex]].dotsOpacity)
      .attr('stroke', d => color(d[groupIndex]))
      .attr('cx', d => x(d[depthIndex]))
      .attr('cy', d => yScale(d[index]));
  }
  const dotsUpdate = chart.selectAll('.circle').data(points);
  dotsUpdate.exit().transition().remove();
  const dotsEnter = dotsUpdate.enter().append('circle')
    .attr('r', 4);
  dotsUpdate.call(plotDots, medianIndex, y);
  dotsEnter.call(plotDots, medianIndex, y);

  const subDotsUpdate = subChart.selectAll('.circle').data(points);
  subDotsUpdate.exit().transition().remove();
  const subDotsEnter = subDotsUpdate.enter().append('circle')
    .attr('r', 4);
  subDotsUpdate.call(plotDots, countIndex, subY);
  subDotsEnter.call(plotDots, countIndex, subY);

  legendBox.attr('viewBox', `0 0 200 ${ly + 10}`); // CHANGE THIS
  // LINES
  function valueline(yIndex, yScale) {
    return line().x(d => x(d[depthIndex]))
      .y(d => yScale(d[yIndex]));
  }
  const datum = nest()
    .key(d => d[groupIndex])
    .entries(points);
  function plotLines(selection, yScale, yIndex) {
    selection.attr('class', d => `line ${d.key}`)
      .attr('stroke', d => color(d.key))
      .attr('opacity', d => curData[d.key].lineOpacity)
      .attr('fill', 'none')
      .attr('d', d => valueline(yIndex, yScale)(d.values));
  }
  const linesUpdate = chart.selectAll('.line').data(datum);
  linesUpdate.exit().transition().remove();
  linesUpdate.enter().append('path').call(plotLines, y, medianIndex);
  linesUpdate.call(plotLines, y, medianIndex);

  const subLinesUpdate = subChart.selectAll('.line').data(datum);
  subLinesUpdate.exit().transition().remove();
  subLinesUpdate.enter().append('path').call(plotLines, subY, countIndex);
  subLinesUpdate.call(plotLines, subY, countIndex);
}

// re-render chart edges, exis, formatting, etc. when selection changes
export default function render(svg, data, category, legend, legendTitle) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 80, right: 50, bottom: 50 };
  const chart = svg.select('#chart');
  const subChart = svg.select('#subChart');

  const { xAxisLabel, yAxisLabel, minX, maxX, minY, maxY, minSubY, maxSubY } = data;

  const xAxis = axisBottom();
  const yAxisChart = axisLeft();
  const yAxisSubChart = axisLeft();

  function setPad(min, max, axis) {
    let pad = (max - min) * 0.03;
    if (Number.isInteger(min) && Number.isInteger(max)) {
      pad = Math.max(Math.round(pad), 1);
      const between = Math.max(3, (max - min) + (2 * pad));
      axis.ticks(Math.min(between, 12), 'd');
    }
    return pad;
  }

  const xPad = setPad(minX, maxX, xAxis);
  const subYPad = setPad(minSubY, maxSubY, yAxisSubChart);

  const x = scaleLinear().domain([minX - xPad, maxX + xPad]).range([0, width]).nice();
  const y = scaleLinear().domain([minY, maxY]).range([height, 0]).nice();
  console.log('subY values: ', subYPad, minSubY, maxSubY);
  const subY = scaleLinear().domain([minSubY - subYPad, maxSubY + subYPad])
    .range([height, 0]).nice();

  xAxis.scale(x);
  yAxisChart.scale(y);
  yAxisSubChart.scale(subY);

  setupXLabel(svg, width, height, xAxisLabel, xAxis);
  const maxLabelY = setupYLabels(svg, height, yAxisLabel, yAxisChart, yAxisSubChart);
  const moveX = Math.max(margin.left, maxLabelY);

  svg.attr('width', width + moveX + margin.right)
    .attr('height', 2 * (height + margin.bottom + margin.top));
  select(svg.node().parentNode).style('width', `${width + moveX + margin.right}px`)
    .style('height', `${2 * (height + margin.bottom + margin.top)}px`);
  chart.attr('transform', `translate(${moveX},${margin.top})`);
  subChart.attr('transform', `translate(${moveX},${height + margin.bottom + margin.top})`);
  renderPlot(svg, data, x, y, subY, category, legend, legendTitle);
}
