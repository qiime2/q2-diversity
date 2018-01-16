/* global columns */
/* global metrics */
import { select } from 'd3';
import state from './state';
import { addMetricPicker, addColumnPicker } from './toolbar';

export default function init() {
  const metric = metrics[0];
  const column = columns[0];
  const body = select('#main');
  const controlsRow = body.select('.controls');
  const plotSvg = body.select('.plotSvg');
  const legendBox = body.select('.legendBoxSvg');
  const legendTitle = body.select('.legendTitle');
  state.initialize(metric, column, controlsRow, plotSvg, legendBox, legendTitle);
  addMetricPicker(controlsRow, metrics, metric);
  if (columns.length > 0) {
    addColumnPicker(controlsRow, columns, column);
  }
}
