import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory20,
  select,
  line,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';
import appendLegendKey from './legend';
import { curData } from './data';

// re-render chart and legend whenever selection changes
function renderPlot(svg, data, x, y, category, legend, legendTitle) {
  const chart = svg.select('g');

  const depthIndex = data.data.columns.indexOf('depth');
  const medianIndex = data.data.columns.indexOf('median');
  let groupIndex = data.data.columns.indexOf('sample-id');
  if (groupIndex === -1) {
    groupIndex = data.data.columns.indexOf(category);
  }
  const points = [data.data.data.sort((a, b) => a[depthIndex] - b[depthIndex])][0];
  const setGroups = new Set(Array.from(points, d => d[groupIndex]));
  const color = scaleOrdinal(schemeCategory20)
    .domain(setGroups);
  const arrGroups = Array.from(setGroups);

  const legendBox = select(legend.node().parentNode);

  legend.selectAll('.legend').remove();
  legendTitle.selectAll('.legend').remove();
  chart.selectAll('.circle').remove();
  chart.selectAll('.line').remove();

  legend.attr('height', arrGroups.length * 20);
  let ly = 0;

  const valueline = line()
    .x(d => x(d[depthIndex]))
    .y(d => y(d[medianIndex]));

  appendLegendKey(legendTitle, 0, 'Select%20All', 10, 'black',
                  color, arrGroups, chart);
  for (const [i, entry] of arrGroups.entries()) {
    ly = (i + 0.5) * 20;
    const subset = points.filter(d => d[groupIndex] === entry)
                    .sort((a, b) => a[depthIndex] - b[depthIndex]);
    curData.appendSeries(entry, subset);
    curData.toggle(entry, 1, 0);
    appendLegendKey(legend, i + 1, entry, ly, color(entry), color,
                    arrGroups, chart);
    const curColor = color(subset);
    chart.append('path')
        .attr('d', valueline(curData.getSeries(entry)))
        .style('stroke', curColor)
        .style('opacity', curData.line(entry))
        .style('fill', 'none')
        .attr('class', 'line')
        .attr('id', `idline${entry}`);
    chart.selectAll('dot')
        .data(curData.getSeries(entry))
      .enter()
        .append('circle')
          .attr('cx', d => x(d[depthIndex]))
          .attr('cy', d => y(d[medianIndex]))
          .attr('r', 4)
          .style('stroke', curColor)
          .style('opacity', curData.dots(entry))
          .style('fill', curColor)
          .attr('class', 'circle')
          .attr('id', d => `idcircle${d[groupIndex]}`);
  }
  legendBox.attr('viewBox', `0 0 200 ${ly + 10}`);
}

// re-render chart edges, exis, formatting, etc. when selection changes
export default function render(svg, data, category, legend, legendTitle) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 70, right: 50, bottom: 50 };
  const chart = svg.select('g');

  const { xAxisLabel, yAxisLabel, minX, maxX, minY, maxY } = data;

  const xAxis = axisBottom();
  const yAxis = axisLeft();

  let pad = (maxX - minX) * 0.03;
  if (Number.isInteger(minX) && Number.isInteger(maxX)) {
    pad = Math.max(Math.round(pad), 1);
    const between = Math.max(3, (maxX - minX) + (2 * pad));
    xAxis.ticks(Math.min(between, 12), 'd');
  }

  const x = scaleLinear().domain([minX - pad, maxX + pad]).range([0, width]).nice();
  const y = scaleLinear().domain([minY, maxY]).range([height, 0]).nice();

  xAxis.scale(x);
  yAxis.scale(y);

  chart.attr('transform', `translate(${margin.left},${margin.top})`);

  setupXLabel(svg, width, height, xAxisLabel, xAxis);
  setupYLabel(svg, height, yAxisLabel, yAxis);

  renderPlot(svg, data, x, y, category, legend, legendTitle);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}
