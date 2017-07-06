import {
  select,
  scaleLinear,
  axisBottom,
  axisLeft,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';
import plotScatter from './scatter';


export function render(svg, data) {
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

  plotScatter(chart, data, x, y);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}

export function stats(body, data) {
  const { stats: { method, testStat, pVal, sampleSize } } = data;
  select('#method').text(method);
  select('#test-stat').text(testStat);
  select('#p-val').html(pVal);
  select('#sample-size').html(sampleSize);
}

export function warnings(body, data) {
  if (data.filtered === null) {
    select('#filtered-samples').style('display', 'none').html(null);
    return;
  }

  const { filtered: { initial, filtered, method } } = data;
  if (initial !== filtered) {
    select('#filtered-samples')
      .style('display', null)
      .html(`Some samples were filtered from the input alpha diversity data
        because they were missing metadata values.<strong>The input alpha
        diversity data contained ${initial} samples, but ${method} correlation
        was computed on only ${filtered} samples.</strong>`);
  } else {
    select('#filtered-samples').style('display', 'none').html(null);
  }
}
